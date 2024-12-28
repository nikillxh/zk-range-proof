use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::RistrettoPoint;
use rand_chacha::ChaCha20Rng;
use rand::SeedableRng;
use rand_core::OsRng;
use sha2::{Sha512, Digest};

pub fn gen_basis_vectors(count: usize, seed: &[u8], domain: &str) -> Vec<RistrettoPoint> {
    let mut rng = ChaCha20Rng::from_seed(
        Sha512::digest(&[seed, domain.as_bytes()].concat()).as_slice()[..32]
            .try_into()
            .expect("Hash output must fit in 32 bytes"),
    );

    println!("Basis vectors generated");

    (0..count)
        .map(|_| RistrettoPoint::random(&mut rng))
        .collect()
}

pub fn gen_scalars(count: usize) -> Vec<Scalar> {
    let mut rng = OsRng; // Secure RNG
    (0..count).map(|_| Scalar::random(&mut rng)).collect()
}

pub fn gen_v_points(mut rng: OsRng) -> Vec<RistrettoPoint> {
    let g = RistrettoPoint::random(&mut rng);
    let b = RistrettoPoint::random(&mut rng);
    println!("Generated G, B");
    vec![g, b]
}
