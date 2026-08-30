use termpdb::render::{
    DEFAULT_CELL_PIXEL_HEIGHT, DEFAULT_CELL_PIXEL_WIDTH, GraphicsBackend, encode_kitty_delete,
    encode_kitty_graphics_rgba, get_terminal_cell_size,
};

#[test]
fn test_encode_kitty_graphics_single_chunk() {
    let rgba = vec![255, 0, 0, 255, 0, 255, 0, 255]; // 2 pixels
    let seq = encode_kitty_graphics_rgba(2, 1, 10, 5, 0, 1, -1, 1, &rgba);

    assert!(seq.starts_with("\x1b[2;1H\x1b_G"));
    assert!(seq.contains("a=T"));
    assert!(seq.contains("f=32"));
    assert!(seq.contains("s=2,v=1"));
    assert!(seq.contains("c=10,r=5"));
    assert!(seq.contains("z=-1"));
    assert!(seq.contains("i=1"));
    assert!(seq.contains("q=2"));
    assert!(seq.contains("m=0"));
    assert!(seq.ends_with("\x1b\\"));
}

#[test]
fn test_encode_kitty_graphics_multi_chunk() {
    // 4000 pixels = 16000 bytes RGBA -> ~21336 bytes Base64 (multiple 4096-byte chunks)
    let rgba = vec![128u8; 4000 * 4];
    let seq = encode_kitty_graphics_rgba(200, 20, 80, 24, 5, 10, -1, 42, &rgba);

    assert!(seq.starts_with("\x1b[11;6H\x1b_Ga=T,f=32,s=200,v=20,c=80,r=24,z=-1,i=42,q=2,m=1;"));
    assert!(seq.contains("\x1b\\"));
    // Intermediate and ending chunk escape sequences
    assert!(seq.contains("\x1b_Gm=1;"));
    assert!(seq.contains("\x1b_Gm=0;"));
    assert!(seq.ends_with("\x1b\\"));
}

#[test]
fn test_encode_kitty_delete() {
    assert_eq!(encode_kitty_delete(Some(1)), "\x1b_Ga=d,d=i,i=1,q=2\x1b\\");
    assert_eq!(encode_kitty_delete(None), "\x1b_Ga=d,d=a,q=2\x1b\\");
    assert_eq!(encode_kitty_delete(Some(2)), "\x1b_Ga=d,d=a,q=2\x1b\\");
}

#[test]
fn test_graphics_backend_toggle_and_default() {
    let mut backend = GraphicsBackend::default();
    assert_eq!(backend, GraphicsBackend::HalfBlock);
    assert!(!backend.is_kitty());

    backend.toggle();
    assert_eq!(backend, GraphicsBackend::Kitty);
    assert!(backend.is_kitty());

    backend.toggle();
    assert_eq!(backend, GraphicsBackend::HalfBlock);
    assert!(!backend.is_kitty());
}

#[test]
fn test_get_terminal_cell_size() {
    let (w, h) = get_terminal_cell_size();
    assert!(w >= 1);
    assert!(h >= 1);
    // In headless test environments, it should fallback to defaults
    assert_eq!(DEFAULT_CELL_PIXEL_WIDTH, 16);
    assert_eq!(DEFAULT_CELL_PIXEL_HEIGHT, 32);
}
