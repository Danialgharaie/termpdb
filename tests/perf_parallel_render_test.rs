use termpdb::render::buffer::Framebuffer;
use termpdb::render::lighting::Lighting;
use termpdb::render::postprocess::{PostProcessConfig, apply_postprocessing};
use termpdb::render::representations::vdw::render_vdw;
use termpdb::render::representations::{LodLevel, RenderContext, build_render_cache};
use termpdb::render::{Camera, ColorScheme, LodMode, Visibility};

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
fn test_parallel_postprocessing_equivalence() {
    let mut fb1 = Framebuffer::new(40, 40);
    let mut fb2 = Framebuffer::new(40, 40);

    for y in 0..40 {
        for x in 0..40 {
            let depth = if (10..30).contains(&x) && (10..30).contains(&y) {
                5.0
            } else {
                15.0
            };
            fb1.set_pixel(x, y, depth, (180, 180, 180));
            fb2.set_pixel(x, y, depth, (180, 180, 180));
        }
    }

    let config = PostProcessConfig {
        outline: true,
        ssao: true,
        outline_threshold: 0.1,
        ssao_radius: 2,
    };

    apply_postprocessing(&mut fb1, &config);
    apply_postprocessing(&mut fb2, &config);

    assert_eq!(fb1.pixels, fb2.pixels);
    assert_eq!(fb1.depth, fb2.depth);
}

#[test]
fn test_vdw_rendering_with_culling() {
    let structure = termpdb::parser::parse_pdb(SAMPLE_1CRN).unwrap();
    let (colors, visible, com, radius, max_vdw) = build_render_cache(
        &structure,
        ColorScheme::Cpk,
        Visibility::default(),
        LodMode::Full,
    );

    let mut camera = Camera::new();
    camera.fit_structure(com, radius);
    let mats = camera.matrices();
    let lighting = Lighting::default();

    let ctx = RenderContext {
        structure: &structure,
        camera: &camera,
        mats,
        lighting: &lighting,
        visibility: Visibility::default(),
        lod: LodLevel::Full,
        colors: &colors,
        visible: &visible,
        com,
        radius,
        max_vdw,
        ribbon_geometry: None,
    };

    let mut fb = Framebuffer::new(60, 40);
    fb.clear((0, 0, 0));
    render_vdw(&ctx, &mut fb);

    let drawn = fb.pixels.iter().filter(|&&p| p != (0, 0, 0)).count();
    assert!(drawn > 100, "VDW mode should render visible atom spheres");
}
