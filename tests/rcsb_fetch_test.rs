//! Offline tests for the RCSB fetch fallback chain.
//!
//! A throwaway std-only HTTP server (`std::net::TcpListener` bound to
//! 127.0.0.1:0) stands in for files.rcsb.org; its endpoints are injected into
//! [`termpdb::parser::rcsb::fetch_pdb_from`] so the multi-URL fallback and the
//! gzip-decompression paths run hermetically, without touching the network.

use std::io::{Read, Write};
use std::net::TcpListener;

use flate2::Compression;
use flate2::write::GzEncoder;

use termpdb::error::TermPdbError;
use termpdb::parser::rcsb::fetch_pdb_from;

/// One canned response: requests for exactly `path` are answered with
/// `status_line`, `content_type`, and `body`. Unrouted paths get a bare 404.
struct Route {
    path: &'static str,
    status_line: &'static str,
    content_type: &'static str,
    body: Vec<u8>,
}

fn route(path: &'static str, content_type: &'static str, body: Vec<u8>) -> Route {
    Route {
        path,
        status_line: "200 OK",
        content_type,
        body,
    }
}

fn not_found(path: &'static str) -> Route {
    Route {
        path,
        status_line: "404 Not Found",
        content_type: "text/plain",
        body: Vec::new(),
    }
}

/// Binds a listener on 127.0.0.1:0 and spawns a thread answering one request
/// at a time from `routes`. Returns `host:port` (the socket is already bound,
/// so connections queue in the backlog before the thread even runs).
///
/// Responses always carry `Content-Length` and `Connection: close`, which is
/// the minimal HTTP/1.x dialect ureq needs from a hand-rolled server.
fn spawn_http_server(routes: Vec<Route>) -> String {
    let listener = TcpListener::bind(("127.0.0.1", 0)).expect("bind test HTTP server");
    let addr = listener.local_addr().expect("local address");
    std::thread::spawn(move || {
        for stream in listener.incoming() {
            let Ok(mut stream) = stream else { continue };
            let Ok(head) = read_request_head(&mut stream) else {
                continue;
            };
            // "GET /path HTTP/1.1" -> "/path"
            let path = head.split_whitespace().nth(1).unwrap_or("");
            let responded = routes.iter().find(|r| r.path == path);
            let (status_line, content_type, body) = match responded {
                Some(r) => (r.status_line, r.content_type, r.body.as_slice()),
                None => ("404 Not Found", "text/plain", &[][..]),
            };
            let head_out = format!(
                "HTTP/1.1 {status_line}\r\nContent-Type: {content_type}\r\n\
                 Content-Length: {}\r\nConnection: close\r\n\r\n",
                body.len()
            );
            let _ = stream.write_all(head_out.as_bytes());
            let _ = stream.write_all(body);
            let _ = stream.flush();
            // Dropping the stream closes the connection per `Connection: close`.
        }
    });
    format!("127.0.0.1:{}", addr.port())
}

/// Reads from the socket until the `\r\n\r\n` end-of-head marker arrives.
fn read_request_head(stream: &mut std::net::TcpStream) -> std::io::Result<String> {
    let mut buf = Vec::new();
    let mut chunk = [0u8; 1024];
    loop {
        let n = stream.read(&mut chunk)?;
        if n == 0 {
            break;
        }
        buf.extend_from_slice(&chunk[..n]);
        if buf.windows(4).any(|w| w == b"\r\n\r\n") {
            break;
        }
    }
    Ok(String::from_utf8_lossy(&buf).into_owned())
}

/// Gzip-compresses `data` exactly as a structure file would be served.
fn gzipped(data: &[u8]) -> Vec<u8> {
    let mut encoder = GzEncoder::new(Vec::new(), Compression::fast());
    encoder.write_all(data).expect("gzip write");
    encoder.finish().expect("gzip finish")
}

/// Minimal valid PDB payload used as the download body.
const SAMPLE_PDB: &str =
    "ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N\nEND\n";

#[test]
fn fetch_falls_back_to_second_url_and_decompresses_gzip() {
    // First endpoint answers 404; second serves gzip-compressed PDB bytes
    // with no Content-Encoding header, forcing fetch_pdb_from's own
    // magic-byte detection + GzDecoder path to do the decompression.
    let base = spawn_http_server(vec![
        not_found("/download/1CRN.pdb.gz"),
        route(
            "/download/1CRN.cif.gz",
            "application/gzip",
            gzipped(SAMPLE_PDB.as_bytes()),
        ),
    ]);
    let urls = vec![
        format!("http://{base}/download/1CRN.pdb.gz"),
        format!("http://{base}/download/1CRN.cif.gz"),
    ];

    let text = fetch_pdb_from("1crn", &urls).expect("second URL must satisfy the fetch");
    assert_eq!(
        text, SAMPLE_PDB,
        "payload must decompress to the original bytes"
    );
}

#[test]
fn fetch_all_urls_404_is_network_error() {
    let base = spawn_http_server(vec![
        not_found("/download/1CRN.pdb.gz"),
        not_found("/download/1CRN.cif.gz"),
        not_found("/mirror/1CRN.pdb.gz"),
    ]);
    let urls = vec![
        format!("http://{base}/download/1CRN.pdb.gz"),
        format!("http://{base}/download/1CRN.cif.gz"),
        format!("http://{base}/mirror/1CRN.pdb.gz"),
    ];

    let err = fetch_pdb_from("1CRN", &urls).expect_err("every URL 404s, so the fetch must fail");
    assert!(
        matches!(err, TermPdbError::NetworkError(_)),
        "expected NetworkError, got: {err:?}"
    );
}

#[test]
fn fetch_connection_refused_is_network_error() {
    // Bind a port, then immediately release it: nothing is listening there.
    let listener = TcpListener::bind(("127.0.0.1", 0)).expect("bind");
    let port = listener.local_addr().expect("addr").port();
    drop(listener);

    let urls = vec![format!("http://127.0.0.1:{port}/download/1CRN.pdb.gz")];
    let err =
        fetch_pdb_from("1CRN", &urls).expect_err("refused connection must surface as an error");
    assert!(
        matches!(err, TermPdbError::NetworkError(_)),
        "expected NetworkError, got: {err:?}"
    );
}

#[test]
fn fetch_empty_url_list_is_network_error() {
    let err = fetch_pdb_from("1CRN", &[]).expect_err("no candidate URLs cannot succeed");
    assert!(
        matches!(err, TermPdbError::NetworkError(_)),
        "expected NetworkError, got: {err:?}"
    );
}

#[test]
fn fetch_rejects_invalid_ids_before_any_request() {
    // The ID is validated up front, so even a reachable-looking URL list is
    // never touched for IDs that could corrupt path semantics.
    for bad in ["", "   ", "1 CRN", "../etc/passwd", "crn;drop"] {
        let urls = vec!["http://127.0.0.1:9/download/x".to_string()];
        let err = fetch_pdb_from(bad, &urls)
            .expect_err("invalid IDs must be rejected before any request");
        assert!(
            matches!(err, TermPdbError::InvalidStructure(_)),
            "expected InvalidStructure for '{bad}', got: {err:?}"
        );
    }
}
