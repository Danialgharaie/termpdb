use termpdb::render::buffer::Framebuffer;

#[test]
fn test_framebuffer_bands_partitioning() {
    let mut fb = Framebuffer::new(80, 48);
    fb.clear((10, 20, 30));

    // Divide into 4 bands of 12 rows each
    let band_height = 12;
    let bands = fb.par_bands_mut(band_height);
    assert_eq!(bands.len(), 4);

    for (i, band) in bands.into_iter().enumerate() {
        assert_eq!(band.width, 80);
        assert_eq!(band.height, 12);
        assert_eq!(band.y_offset, i * 12);
        assert_eq!(band.pixels.len(), 80 * 12);
        assert_eq!(band.depth.len(), 80 * 12);
    }
}

#[test]
fn test_framebuffer_band_pixel_writes() {
    let mut fb = Framebuffer::new(40, 20);
    fb.clear((0, 0, 0));

    let bands = fb.par_bands_mut(10);
    assert_eq!(bands.len(), 2);

    for mut band in bands {
        // Write a test pixel in local coordinates (5, 5) of the band
        band.set_pixel(5, 5, 2.5, (255, 128, 64));
    }

    // Band 0 local (5, 5) -> global (5, 5)
    assert_eq!(fb.get_pixel(5, 5), Some((255, 128, 64)));
    assert_eq!(fb.get_depth(5, 5), Some(2.5));

    // Band 1 local (5, 5) -> global (5, 15)
    assert_eq!(fb.get_pixel(5, 15), Some((255, 128, 64)));
    assert_eq!(fb.get_depth(5, 15), Some(2.5));
}
