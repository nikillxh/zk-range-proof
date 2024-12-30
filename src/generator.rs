use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::RistrettoPoint;
use rand_chacha::ChaCha20Rng;
use rand::SeedableRng;
use rand_core::OsRng;
use sha2::{Sha512, Digest};
use crate::operations::scalarize;

pub struct GlobalPoints {
    g_basis: Vec<RistrettoPoint>,
    h_basis: Vec<RistrettoPoint>,
    g_i: RistrettoPoint,
    b_i: RistrettoPoint,
}

impl GlobalPoints {
    pub fn gen_global(range: usize) -> Self {
        let seed = b"G and H basis seed";
        let g_basis = gen_basis_vectors(range, seed, "g_basis");
        let h_basis = gen_basis_vectors(range, seed, "h_basis");
    
        let mut rng = OsRng;
        let g_i = RistrettoPoint::random(&mut rng);
        let b_i = RistrettoPoint::random(&mut rng);
        println!("Generated G_basis, H_basis, G, B");

        Self {
            g_basis,
            h_basis,
            g_i,
            b_i,
        }
    }

    pub fn g_basis(&self) -> Vec<RistrettoPoint> {
        self.g_basis.clone()
    }

    pub fn h_basis(&self) -> Vec<RistrettoPoint> {
        self.h_basis.clone()
    }

    pub fn g_i(&self) -> RistrettoPoint {
        self.g_i.clone()
    }

    pub fn b_i(&self) -> RistrettoPoint {
        self.b_i.clone()
    }
}

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

pub fn gen_scalars(count: usize, y: u64) -> Vec<Scalar> {
    let mut rng = ChaCha20Rng::seed_from_u64(y);
    (0..count).map(|_| Scalar::random(&mut rng)).collect()
}

pub fn n2_gen(range: usize) -> Vec<Scalar> {
    scalarize(&mut vec![2i64; range].iter().enumerate().map(|(x,i)| i.pow(x as u32)).collect())
}
