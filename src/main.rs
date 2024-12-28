use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::RistrettoPoint;
use rand::rngs::OsRng;

mod prover;
mod generator;
mod operations;

fn main() {
    let mut rng = OsRng;

    // Generators for the Pedersen commitment
    let g = RistrettoPoint::random(&mut rng);
    let h = RistrettoPoint::random(&mut rng);

    println!("g: {:?}", g);
    println!("h: {:?}", h);

    // Secret value and blinding factor
    let value = Scalar::from(42u64);
    let blinding_factor = Scalar::random(&mut rng);

    // Pedersen commitment: C = v*G + r*H
    let commitment = g * value + h * blinding_factor;

    println!("Pedersen Commitment: {:?}", commitment);
}
