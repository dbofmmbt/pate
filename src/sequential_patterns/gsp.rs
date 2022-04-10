use super::{Dataset, SequentialPatternMiner};

pub struct Gsp;

impl SequentialPatternMiner for Gsp {
    fn mine(&mut self, dataset: &Dataset, sup_min: f32) -> Vec<super::SequentialPattern> {
        let frequent_items = frequent_set_level_1(&dataset);

        todo!()
    }
}

fn frequent_set_level_1(dataset: &super::Dataset) -> () {
    todo!()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn it_works() {}
}
