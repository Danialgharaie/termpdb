use termpdb::render::braille::BrailleBuffer;

#[test]
fn test_braille_dot_mapping() {
    // 1 cell width, 1 cell height => 2x4 subpixels
    let mut buf = BrailleBuffer::new(2, 4);

    // Set dot (0, 0) -> bit 0 (0x01) -> '⠁'
    buf.set_subpixel(0, 0, 1.0, (255, 255, 255));
    let (ch1, _) = buf.cell_at(0, 0);
    assert_eq!(ch1, '⠁');

    // Set dot (1, 3) -> bit 7 (0x80) -> 0x2800 | 0x01 | 0x80 = 0x2881 -> '⢁'
    buf.set_subpixel(1, 3, 1.0, (255, 255, 255));
    let (ch2, _) = buf.cell_at(0, 0);
    assert_eq!(ch2, char::from_u32(0x2800 | 0x01 | 0x80).unwrap());
}

#[test]
fn test_braille_line_rasterization() {
    let mut buf = BrailleBuffer::new(20, 40);
    buf.draw_line_3d((0.0, 0.0, 1.0), (19.0, 39.0, 1.0), (100, 200, 255));

    // Top-left cell should have dots
    let (ch_tl, _) = buf.cell_at(0, 0);
    assert_ne!(ch_tl, ' ');
    assert_ne!(ch_tl, '⠀'); // Not blank Braille (0x2800)

    // Bottom-right cell should have dots
    let (ch_br, _) = buf.cell_at(9, 9);
    assert_ne!(ch_br, ' ');
    assert_ne!(ch_br, '⠀');
}
