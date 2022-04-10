mod apriori;
pub use apriori::Apriori;

pub type Dataset = Vec<Transaction>;
pub type Transaction = ();

pub struct ItemSet;

pub trait ItemSetMiner {
    fn mine(&mut self, dataset: &Dataset, sup_min: f32);
}
