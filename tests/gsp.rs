use pate::sequential_patterns::{Gsp, SequentialPatternMiner, Dataset};

#[test]
fn gsp_works() {
    

    
}

fn check(input: impl Into<Dataset>, expected_output: ()) {
    let mut miner = Gsp;
    miner.mine(todo!(), todo!());
}
