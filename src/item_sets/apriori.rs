use super::{Dataset, ItemSetMiner};

pub struct Apriori;

impl ItemSetMiner for Apriori {
    fn mine(&mut self, dataset: &Dataset, sup_min: f32) {
        let mut current_set = frequent_set_level_1(dataset);
        let mut k = 1;

        while !current_set.is_empty() {
            k += 1;
            let candidates = self.generate_candidates(&current_set, k);

            for transaction in dataset.iter() {}
        }
    }
}

fn frequent_set_level_1(dataset: &Dataset) -> Vec<()> {
    todo!()
}
impl Apriori {
    fn generate_candidates(&self, current_set: &[()], k: i32) -> () {
        todo!()
    }
}
