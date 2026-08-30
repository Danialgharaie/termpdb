use termpdb::render::{Framebuffer, PixelColor};

#[test]
fn test_framebuffer_to_rgba_bytes() {
    let mut fb = Framebuffer::new(2, 2);
    fb.set_pixel(0, 0, 1.0, (255, 0, 0));
    fb.set_pixel(1, 0, 1.0, (0, 255, 0));
    fb.set_pixel(0, 1, 1.0, (0, 0, 255));
    fb.set_pixel(1, 1, 1.0, (255, 255, 255));

    let rgba = fb.to_rgba_bytes();
    assert_eq!(rgba.len(), 2 * 2 * 4);
    // (0,0): R=255, G=0, B=0, A=255
    assert_eq!(&rgba[0..4], &[255, 0, 0, 255]);
    // (1,0): R=0, G=255, B=0, A=255
    assert_eq!(&rgba[4..8], &[0, 255, 0, 255]);
    // (0,1): R=0, G=0, B=255, A=255
    assert_eq!(&rgba[8..12], &[0, 0, 255, 255]);
    // (1,1): R=255, G=255, B=255, A=255
    assert_eq!(&rgba[12..16], &[255, 255, 255, 255]);
}

#[test]
fn test_framebuffer_write_rgba_bytes_reuses_allocation() {
    let mut fb = Framebuffer::new(3, 1);
    let c1: PixelColor = (10, 20, 30);
    let c2: PixelColor = (40, 50, 60);
    let c3: PixelColor = (70, 80, 90);
    fb.set_pixel(0, 0, 1.0, c1);
    fb.set_pixel(1, 0, 1.0, c2);
    fb.set_pixel(2, 0, 1.0, c3);

    let mut buf = vec![99u8; 100]; // Pre-existing data
    fb.write_rgba_bytes(&mut buf);

    assert_eq!(buf.len(), 3 * 4);
    assert_eq!(&buf[0..4], &[10, 20, 30, 255]);
    assert_eq!(&buf[4..8], &[40, 50, 60, 255]);
    assert_eq!(&buf[8..12], &[70, 80, 90, 255]);
}

#[test]
fn test_framebuffer_to_rgba_bytes_empty() {
    let fb = Framebuffer::new(0, 0);
    let rgba = fb.to_rgba_bytes();
    assert_eq!(rgba.len(), 0);
}
