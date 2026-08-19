use thiserror::Error;

pub type Result<T> = std::result::Result<T, TermPdbError>;

#[derive(Error, Debug)]
pub enum TermPdbError {
    #[error("I/O error: {0}")]
    Io(#[from] std::io::Error),

    #[error("Parse error: {0}")]
    ParseError(String),

    #[error("Network error: {0}")]
    NetworkError(String),

    #[error("Invalid structure data: {0}")]
    InvalidStructure(String),

    #[error("Render error: {0}")]
    RenderError(String),

    #[error("Math error: {0}")]
    MathError(String),

    #[error("Error: {0}")]
    Other(String),
}
