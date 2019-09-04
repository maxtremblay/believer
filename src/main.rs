use believer::CodeGenerator;

fn main() {
    let generator = CodeGenerator::new(1.0, 0.9, 1.0, 4, 20, 15);
    let s = 2;
    let matrix = generator.generate();
    matrix.rows_iter().for_each(|check| {
        println!("{:?}", check.positions());
    })
}
