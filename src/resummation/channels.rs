use crate::models::hubbard::HubbardModel;

/// Particle-particle ladder resummation channel
#[derive(Debug, Clone)]
pub struct PPLadder {
    pub model: HubbardModel,
}

impl PPLadder {
    pub fn new(model: HubbardModel) -> Self {
        Self { model }
    }
}
