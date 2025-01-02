use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::RistrettoPoint;
use rand::Rng;
use rand_core::OsRng;

use crate::bullerproof::fold_vector;
use crate::generator::{n2_gen, gen_scalars, GlobalPoints};
use crate::operations::{inner_product, inv_vector, points_hadamard_multiply, scalarize, vec_scalar_mul, vector_add};
use crate::prover::{Polycommitment, T1T2commitment};

pub struct Generatives {
    y: u64,
    z: u64,
    u: Scalar,
    yn: Vec<Scalar>,
    y_inv_h: Vec<RistrettoPoint>,
}

pub struct LinearVerify {
    lu: Vec<Scalar>,
    ru: Vec<Scalar>,
    tu: Scalar,
    eqn2lhs: RistrettoPoint,
    eqn2rhs: RistrettoPoint,
    eqn3lhs: RistrettoPoint,
    eqn3rhs: RistrettoPoint,
}

pub struct BulletVerify {
    commit_a: RistrettoPoint,
    commit_s: RistrettoPoint,
    commit_v: RistrettoPoint,
    commit_c: RistrettoPoint,
    commit_p: RistrettoPoint,
    commit_t1: RistrettoPoint,
    commit_t2: RistrettoPoint,
    g_basis_fold: Vec<RistrettoPoint>,
    h_basis_fold: Vec<RistrettoPoint>,
    g_i: RistrettoPoint,
    u_random: Scalar,
    pi_lr: Scalar,
    pi_t: Scalar,
    tu: Scalar,
    z: Scalar,
}

impl BulletVerify {
    pub fn init([left, right]: [RistrettoPoint; 2], asv: [RistrettoPoint; 3], [commit_t1, commit_t2]: [RistrettoPoint; 2], data: (RistrettoPoint, [Scalar; 3]), points: GlobalPoints, y_inv_h: Vec<RistrettoPoint>, z: u64, count: usize) -> Self {
        let u_random = Scalar::random(&mut OsRng);
        let commit_c = data.0;
        let [tu, pi_lr, pi_t] = data.1;
        let mut commit_p = commit_c + (tu * points.g_i());
        commit_p = (left * u_random * u_random) + (right * u_random.invert() * u_random.invert()) + commit_p;
        let g_basis_fold = fold_vector(points.g_basis().clone(), u_random.invert());
        let h_basis_fold = fold_vector(y_inv_h, u_random);
        let z = Scalar::from(z);
        Self {
            commit_a: asv[0],
            commit_s: asv[1],
            commit_v: asv[2],
            commit_c,
            commit_p,
            commit_t1,
            commit_t2,
            g_basis_fold,
            h_basis_fold,
            g_i: points.g_i(),
            u_random,
            pi_lr,
            pi_t,
            tu,
            z
        }
    }

    pub fn compute(&mut self, [left, right]: [RistrettoPoint; 2]) -> () {
        self.commit_p = (left * self.u_random * self.u_random) + (right * self.u_random.invert() * self.u_random.invert()) + self.commit_p;
        self.g_basis_fold = fold_vector(self.g_basis_fold.clone(), self.u_random);
        self.h_basis_fold = fold_vector(self.h_basis_fold.clone(), self.u_random);
    }

    pub fn u_update(&mut self) -> Scalar {
        self.u_random = Scalar::random(&mut OsRng);
        self.u_random
    }

    pub fn verify(&mut self, [a, b]: [Vec<Scalar>; 2], count: usize, points: GlobalPoints, gen: Generatives) -> () {
        let n2 = n2_gen(count);
        let delta: Scalar = ((Scalar::from(gen.z()) - Scalar::from(gen.z().pow(2))) * gen.yn().iter().sum::<Scalar>()) 
            - (Scalar::from(gen.z().pow(3)) * n2.iter().sum::<Scalar>());

        let eqn1lhs = self.commit_p;
        let eqn1rhs = (a[0] * self.g_basis_fold[0]) + (b[0] * self.h_basis_fold[0]) + (a[0] * b[0] * self.g_i);
        let eqn2lhs = self.commit_a + (self.commit_s * self.u_random) + inner_product(&vec![Scalar::from(0u64) - self.z; count], &points.g_basis())
            + inner_product(&vector_add(&vec_scalar_mul(&gen.yn(), &Scalar::from(self.z)), &vec_scalar_mul(&n2, &Scalar::from(gen.z()))), 
            &gen.y_inv_h());
        let eqn2rhs = self.commit_c;
        let eqn3lhs = (self.tu * points.g_i()) + (self.pi_t * points.b_i());
        let eqn3rhs = (self.commit_v * Scalar::from(gen.z().pow(2))) + (delta * points.g_i()) + (self.commit_t1 * gen.u()) + (self.commit_t2 * gen.u());

        assert_eq!(eqn1lhs, eqn1rhs, "Final, P != aG + bH + ab(G_i)");
        assert_eq!(eqn2lhs, eqn2rhs, "Final, A, S verification failed");
        assert_eq!(eqn3lhs, eqn3rhs, "Final, V, T1, T2 verification failed");
    }
}

impl LinearVerify {
    pub fn init_linear(count: usize, points: GlobalPoints, prover: Polycommitment, asv: [RistrettoPoint; 3], gen: Generatives, t_commit: T1T2commitment) -> Self {
        let [commit_a, commit_s, commit_v] = asv;
        let n2 = n2_gen(count);
        let delta: Scalar = ((Scalar::from(gen.z()) - Scalar::from(gen.z().pow(2))) * gen.yn().iter().sum::<Scalar>()) 
            - (Scalar::from(gen.z().pow(3)) * n2.iter().sum::<Scalar>());

        let eqn2lhs = commit_a + (commit_s * gen.u()) + inner_product(&scalarize(&mut vec![-1 * gen.z() as i64; count]), &points.g_basis())
            + inner_product(&vector_add(&vec_scalar_mul(&gen.yn(), &Scalar::from(gen.z())), &vec_scalar_mul(&n2, &Scalar::from(gen.z()))), 
            &gen.y_inv_h());
        let eqn2rhs = inner_product(&prover.lu(), &points.g_basis()) + inner_product(&prover.ru(), &gen.y_inv_h()) 
            + (prover.pi_lr() * points.b_i());
        
        let eqn3lhs = (prover.tu() * points.g_i()) + (prover.pi_t() * points.b_i());
        let eqn3rhs = (commit_v * Scalar::from(gen.z().pow(2))) + (delta * points.g_i()) + (t_commit.commit_t1() * gen.u()) + (t_commit.commit_t2() * gen.u());

        println!("Verification setup initialized");

        Self {
            lu: prover.lu(),
            ru: prover.ru(),
            tu: prover.tu(),
            eqn2lhs,
            eqn2rhs,
            eqn3lhs,
            eqn3rhs,
        }
    }

    pub fn gen_u() -> Scalar {
        Scalar::random(&mut OsRng)
    }

    pub fn verify() -> () {

    }
}

impl Generatives {
    pub fn init(count: usize, points: &GlobalPoints) -> Self {
        let y = rand::thread_rng().gen();
        let z = rand::thread_rng().gen();
        let u = Scalar::random(&mut OsRng);
        let yn = gen_scalars(count, y);
        let yinv = inv_vector(&yn);
        let y_inv_h = points_hadamard_multiply(&yinv, &points.h_basis());
        println!("Generated y, z, y_inv_H");

        Self {
            y,
            z,
            u,
            yn,
            y_inv_h,
        }
    }

    pub fn to_prover_yz(&self) -> [u64; 2] {
        [self.y, self.z]
    }

    pub fn u(&self) -> Scalar {
        self.u
    }

    pub fn z(&self) -> u64 {
        self.z
    }

    pub fn y(&self) -> u64 {
        self.y
    }

    pub fn yn(&self) -> Vec<Scalar> {
        self.yn.clone()
    } 

    pub fn y_inv_h(&self) -> Vec<RistrettoPoint> {
        self.y_inv_h.clone()
    }
}

