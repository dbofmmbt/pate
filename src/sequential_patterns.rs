mod gsp;
pub use gsp::Gsp;

pub struct Dataset;

pub struct SequentialPattern;

pub trait SequentialPatternMiner {
    fn mine(&mut self, dataset: &Dataset, sup_min: f32) -> Vec<SequentialPattern>;
}
