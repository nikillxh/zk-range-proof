use bullerproof::prove_commitments_log;
use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::RistrettoPoint;
use operations::{to_bin, vector_sub, scalarize};
use generator::GlobalPoints;
use prover::{ASVcommitment, BulletProof, Polycommitment, Salts, T1T2commitment};
use rand::rngs::OsRng;
use verifier::{BulletVerify, Generatives, LinearVerify};

mod prover;
mod generator;
mod operations;
mod verifier;
mod bullerproof;

fn main() {
    // Prover
    let prover_value: u64 = 63;
    let al = scalarize(&mut to_bin(prover_value));
    let range = 7;
    let ar = vector_sub( &al, &vec![Scalar::from(1u32)]);
    let salt = Salts::init();

    // Global
    let points= GlobalPoints::gen_global(range);

    // Prover
    let asv = ASVcommitment::compute(prover_value, range, &salt, &points);

    // Verifier
    let gen = Generatives::init(range, &points);
    let [y, z] = [gen.y(), gen.z()];

    // Prover
    let t1t2 = T1T2commitment::init(&salt, &asv, y, z, range, &points);
    let [commit_t1, commit_t2] = [t1t2.commit_t1(), t1t2.commit_t2()];

    // Verifier
    let u = gen.u();

    // Prover
    let poly = Polycommitment::compute(u, salt, &asv, &t1t2, y, z, range, &points);


    // Linear Verification

    // Verifier
    let linear_verifier = LinearVerify::init_linear(range, &points, &poly, asv.to_verifier(), &gen, &t1t2);

    linear_verifier.verify();


    // Bulletproof Verification

    // Prover
    let [left, right] = BulletProof::compute_diagonal([poly.lu(), poly.ru()], [&mut points.g_basis(), &mut poly.y_inv_h()], points.g_i());

    // Verifier
    let mut verifier = BulletVerify::init([left, right], asv.to_verifier(), [commit_t1, commit_t2], poly.bullet_verifier(), &points, gen.y_inv_h(), z, range);

    // Prover
    let mut prover = BulletProof::init(verifier.u_random(), [left, right], &poly, &points);

    // Prove commitments log
    prove_commitments_log(range, &points, &mut prover, &mut verifier, &gen);


    
}
