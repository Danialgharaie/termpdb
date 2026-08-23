use termpdb::parser::parse_pdb;
use termpdb::render::buffer::Framebuffer;
use termpdb::render::export::{downsample_rgba, render_supersampled, render_svg, write_png};
use termpdb::render::{ColorScheme, LodMode, RenderMode, Visibility};

const SAMPLE_1CRN: &str = r#"HEADER    PLANT SEED PROTEIN                     30-APR-81   1CRN
ATOM      1  N   THR A   1      17.047  14.099   3.625  1.00 13.79           N
ATOM      2  CA  THR A   1      16.967  12.784   4.338  1.00 10.80           C
ATOM      3  C   THR A   1      15.685  12.755   5.133  1.00  9.19           C
ATOM      4  O   THR A   1      15.268  13.825   5.594  1.00  9.85           O
ATOM      5  CB  THR A   1      18.170  12.703   5.337  1.00 13.02           C
ATOM      6  OG1 THR A   1      19.334  12.829   4.463  1.00 15.06           O
ATOM      7  CG2 THR A   1      18.150  11.546   6.304  1.00 19.25           C
ATOM      8  N   THR A   2      15.115  11.555   5.265  1.00  7.81           N
ATOM      9  CA  THR A   2      13.856  11.469   6.066  1.00  7.14           C
ATOM     10  C   THR A   2      14.164  10.785   7.379  1.00  5.84           C
ATOM     11  O   THR A   2      14.993   9.862   7.448  1.00  6.94           O
TER      12      THR A   2
END
"#;

#[test]
fn test_downsample_rgba_box_filter() {
    // High-res framebuffer: 4x4, downsampled to 2x2 with SSAA = 2
    let mut fb = Framebuffer::new(4, 4);
    fb.clear((0, 0, 0));

    // Top-left 2x2 block in high-res:
    // set 3 out of 4 pixels to red (200, 0, 0)
    fb.set_pixel(0, 0, 5.0, (200, 0, 0));
    fb.set_pixel(1, 0, 5.0, (200, 0, 0));
    fb.set_pixel(0, 1, 5.0, (200, 0, 0));
    // pixel (1, 1) remains background (depth is infinite)

    let rgba = downsample_rgba(&fb, 2, 2, 2);
    assert_eq!(rgba.len(), 2 * 2 * 4);

    // Pixel (0, 0) in 2x2 output:
    // Color average of 3 pixels = (200, 0, 0)
    // Alpha = (3 * 255) / 4 = 191
    assert_eq!(rgba[0], 200);
    assert_eq!(rgba[1], 0);
    assert_eq!(rgba[2], 0);
    assert_eq!(rgba[3], 191);

    // Pixel (1, 1) in 2x2 output (unpainted area):
    // Alpha should be 0
    let o = (2 + 1) * 4;
    assert_eq!(rgba[o + 3], 0);
}

#[test]
fn test_supersampled_render_modes_and_ssaa_factors() {
    let mut structure = parse_pdb(SAMPLE_1CRN).expect("Failed to parse 1CRN");
    structure.build_bonds();

    let width = 64;
    let height = 64;

    for ssaa in [1, 2, 4] {
        for mode in [
            RenderMode::Ribbon,
            RenderMode::BallAndStick,
            RenderMode::Vdw,
            RenderMode::Trace,
        ] {
            let rgba = render_supersampled(
                &structure,
                mode,
                ColorScheme::Cpk,
                width,
                height,
                ssaa,
                Visibility::ALL,
                LodMode::Full,
            );

            assert_eq!(
                rgba.len(),
                width * height * 4,
                "RGBA output size should match width * height * 4 for mode {mode:?} at SSAA {ssaa}"
            );

            // Count drawn pixels (alpha > 0)
            let drawn_count = rgba.chunks_exact(4).filter(|chunk| chunk[3] > 0).count();
            assert!(
                drawn_count > 50,
                "Mode {mode:?} with SSAA {ssaa} should render visible pixels, got {drawn_count}"
            );

            if ssaa > 1 {
                // In SSAA modes, antialiased edge pixels should exist where 0 < alpha < 255
                let antialiased_pixels = rgba
                    .chunks_exact(4)
                    .filter(|chunk| chunk[3] > 0 && chunk[3] < 255)
                    .count();
                assert!(
                    antialiased_pixels > 0,
                    "SSAA {ssaa} should generate smoothed antialiased boundary pixels for {mode:?}"
                );
            }
        }
    }
}

#[test]
fn test_png_write_and_header_validity() {
    let mut structure = parse_pdb(SAMPLE_1CRN).expect("Failed to parse 1CRN");
    structure.build_bonds();

    let width = 48;
    let height = 48;
    let rgba = render_supersampled(
        &structure,
        RenderMode::Ribbon,
        ColorScheme::Rainbow,
        width,
        height,
        2,
        Visibility::ALL,
        LodMode::Full,
    );

    let temp_path = std::env::temp_dir().join("termpdb_test_supersampled.png");
    let temp_str = temp_path.to_str().unwrap();

    write_png(temp_str, &rgba, width as u32, height as u32).expect("write_png should succeed");

    let bytes = std::fs::read(&temp_path).expect("Failed to read back generated PNG");
    assert!(
        bytes.len() > 100,
        "PNG file should contain encoded image data"
    );

    // Standard PNG magic bytes: 0x89 0x50 0x4E 0x47 0x0D 0x0A 0x1A 0x0A
    let png_magic = [0x89, 0x50, 0x4E, 0x47, 0x0D, 0x0A, 0x1A, 0x0A];
    assert_eq!(
        &bytes[0..8],
        &png_magic,
        "File must have valid PNG magic signature"
    );

    let _ = std::fs::remove_file(&temp_path);
}

#[test]
fn test_render_svg_output() {
    let mut structure = parse_pdb(SAMPLE_1CRN).expect("Failed to parse 1CRN");
    structure.build_bonds();

    let svg = render_svg(
        &structure,
        RenderMode::BallAndStick,
        ColorScheme::Cpk,
        200,
        200,
        Visibility::ALL,
        LodMode::Full,
    );

    assert!(svg.starts_with("<svg"), "SVG must start with <svg root tag");
    assert!(
        svg.ends_with("</svg>\n") || svg.ends_with("</svg>"),
        "SVG must end with </svg>"
    );
    assert!(
        svg.contains("<circle") || svg.contains("<line"),
        "SVG must contain rendered primitives"
    );
}
