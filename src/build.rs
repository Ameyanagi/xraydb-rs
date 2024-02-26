use std::fs;
fn main() {
    fs::write("test.txt", "this is a test");
    println!("This is a build time test");
}
