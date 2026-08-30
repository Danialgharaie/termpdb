use termpdb::model::Structure;
use termpdb::render::{ColorScheme, LodMode, RenderMode, Visibility, export_gif};

#[test]
fn test_export_gif_generates_valid_gif_file() {
    let structure = Structure::default();
    let tmp_path = std::env::temp_dir().join("test_export_turntable.gif");
    let res = export_gif(
        &structure,
        RenderMode::Ribbon,
        ColorScheme::Rainbow,
        80,
        60,
        1,
        3,
        10,
        Visibility::default(),
        LodMode::Auto,
        tmp_path.to_str().unwrap(),
    );
    assert!(res.is_ok());
    assert!(tmp_path.exists());
    let bytes = std::fs::read(&tmp_path).unwrap();
    // GIF magic header: GIF89a or GIF87a
    assert!(bytes.starts_with(b"GIF89a") || bytes.starts_with(b"GIF87a"));
    let _ = std::fs::remove_file(tmp_path);
}
