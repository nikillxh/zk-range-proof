use curve25519_dalek::scalar::Scalar;
use curve25519_dalek::ristretto::RistrettoPoint;
use rand::rngs::OsRng;
use rand::Rng;
use crate::bullerproof::{fold_scalar, fold_vector};
use crate::operations::{points_hadamard_multiply, inv_vector, diagonal_ss_sum, diagonal_sv_sum, diagonal_vs_sum, hadamard_multiply, inner_product, scalarize, to_bin, vec_scalar_mul, vector_add, vector_sub};
use crate::generator::{gen_scalars, n2_gen, GlobalPoints};

pub struct ASVcommitment {
    commit_a: RistrettoPoint,
    commit_s: RistrettoPoint,
    commit_v: RistrettoPoint,
    al: Vec<Scalar>,
    ar: Vec<Scalar>,
    sl: Vec<Scalar>,
    sr: Vec<Scalar>,
}

pub struct Salts {
    salt_alpha: Scalar,
    salt_beta: Scalar,
    salt_gamma: Scalar,
    salt_tau1: Scalar,
    salt_tau2: Scalar,
}

pub struct Polycommitment {
    l: Vec<Scalar>,
    r: Vec<Scalar>,
    t: Scalar,
    pi_lr: Scalar,
    pi_t: Scalar,
    commit_c: RistrettoPoint,
    y_inv_h: Vec<RistrettoPoint>
}

pub struct T1T2commitment {
    commit_t1: RistrettoPoint,
    commit_t2: RistrettoPoint,
    tx: [Scalar; 3]
}

pub struct BulletProof {
    left: RistrettoPoint,
    right: RistrettoPoint,
    commit_p: RistrettoPoint,
    g_basis_fold: Vec<RistrettoPoint>,
    h_basis_fold: Vec<RistrettoPoint>,
    g_i: RistrettoPoint,
    a: Vec<Scalar>,
    b: Vec<Scalar>,
    u_verifier: Scalar,
}

impl BulletProof {
    pub fn init(u_random: Scalar, [left, right]: [RistrettoPoint; 2], poly: &Polycommitment, points: &GlobalPoints) -> Self {
        let commit_p = (left * u_random * u_random) + (right * u_random.invert() * u_random.invert()) + poly.commit_c() + (poly.tu() * points.g_i());
        let g_basis_fold = fold_vector(&mut points.g_basis().clone(), u_random.invert());
        let h_basis_fold = fold_vector(&mut poly.y_inv_h(), u_random);
        let a = fold_scalar(&mut poly.lu(), u_random);
        let b = fold_scalar(&mut poly.ru(), u_random.invert());
        Self {
            left,
            right,
            commit_p,
            g_basis_fold,
            h_basis_fold,
            g_i: points.g_i(),
            a,
            b,
            u_verifier: u_random,
        }
    }

    pub fn compute_diagonal([mut left, mut right]: [Vec<Scalar>; 2], [mut g_basis, mut h_basis]: [&mut Vec<RistrettoPoint>; 2], g_i: RistrettoPoint) -> [RistrettoPoint; 2] {
        let new_left = (diagonal_ss_sum(&mut left, &mut right) * g_i) +
            (diagonal_sv_sum(&mut left, &mut g_basis)) + (diagonal_sv_sum(&mut right, &mut h_basis));
        let new_right = (diagonal_ss_sum(&mut right, &mut left) * g_i) +
            (diagonal_vs_sum(&mut g_basis, &mut left)) + (diagonal_vs_sum(&mut h_basis, &mut right));
        
        [new_left, new_right]
    }

    pub fn a_fold(&mut self) -> () {
        self.a = fold_scalar(&mut self.a(), self.u_verifier);
    }

    pub fn b_fold(&mut self) -> () {
        self.b = fold_scalar(&mut self.b(), self.u_verifier.invert());
    }

    pub fn compute(&mut self, u_random: Scalar) -> () {
        self.u_verifier = u_random;
        self.commit_p = (self.left * self.u_verifier * self.u_verifier) + (self.right * self.u_verifier.invert() * self.u_verifier.invert()) + self.commit_p;
        self.g_basis_fold = fold_vector(&mut self.g_basis_fold.clone(), self.u_verifier.invert());
        self.h_basis_fold = fold_vector(&mut self.h_basis_fold.clone(), self.u_verifier);
        self.a_fold();
        self.b_fold();
    }

    pub fn a(&self) -> Vec<Scalar> {
        self.a.clone()
    }

    pub fn b(&self) -> Vec<Scalar> {
        self.b.clone()
    }

    pub fn update_diagonals(&mut self, [left, right]: [RistrettoPoint; 2]) -> () {
        self.left = left;
        self.right = right;
    }

    pub fn h_basis(&self) -> Vec<RistrettoPoint> {
        self.h_basis_fold.clone()
    }

    pub fn g_basis(&self) -> Vec<RistrettoPoint> {
        self.g_basis_fold.clone()
    }

    pub fn commit_p(&self) -> RistrettoPoint {
        self.commit_p
    }
}

impl T1T2commitment {
    pub fn init(salt: &Salts, asv: &ASVcommitment, y: u64, z: u64, count: usize, points: &GlobalPoints) -> Self {
        let yn = gen_scalars(count, y);
        let n2 = n2_gen(count);
        let z2 = Scalar::from(z) * Scalar::from(z);
        let [al, ar, sl, sr] = asv.polynomial_const();
        
        let t0: Scalar= hadamard_multiply(&vector_sub(&al, &vec![Scalar::from(z); count]),
            &vector_add(&vector_add(&hadamard_multiply(&yn, &ar), &vec_scalar_mul(&yn, &Scalar::from(z))), 
            &vec_scalar_mul(&n2, &z2))).iter().sum();
        let t1: Scalar = hadamard_multiply(&vector_sub(&al, &vec![Scalar::from(z); count]), &hadamard_multiply(&yn, &sr)).iter().sum::<Scalar>()
            + hadamard_multiply(&vector_add(&hadamard_multiply(&yn, &vector_add(&ar, &vec![Scalar::from(z); count])),
            &vec_scalar_mul(&n2, &z2)), &sl).iter().sum::<Scalar>();
        let t2: Scalar = hadamard_multiply(&sl, &hadamard_multiply(&yn, &sr)).iter().sum();

        let commit_t1 = (t1 * points.g_i()) + (salt.tau1() * points.b_i());
        let commit_t2 = (t2 * points.g_i()) + (salt.tau2() * points.b_i());
        println!("Generated T1, T2 commitments");

        Self {
            commit_t1,
            commit_t2,
            tx: [t0, t1, t2]
        }
    }

    pub fn commit_t1(&self) -> RistrettoPoint {
        self.commit_t1
    }

    pub fn commit_t2(&self) -> RistrettoPoint {
        self.commit_t2
    }

    pub fn access_tx(&self) -> [Scalar; 3] {
        self.tx
    }
}

impl Polycommitment {
    pub fn compute(u: Scalar, salt: Salts, asv: &ASVcommitment, tx: &T1T2commitment, y: u64, z: u64, count: usize, points: &GlobalPoints) -> Self {
        let yn = gen_scalars(count, y);
        let y_inv_h = points_hadamard_multiply(&inv_vector(&yn), &points.h_basis());
        let n2 = n2_gen(count);
        let z2 = Scalar::from(z) * Scalar::from(z);
        let [al, ar, sl, sr] = asv.polynomial_const();
        let l = vector_add(&vector_sub(&al, &vec![Scalar::from(z); count]), &vec_scalar_mul(&sl, &u));
        let r = vector_add(&vector_add(&hadamard_multiply(&yn, &vector_add(&ar, &vec![Scalar::from(z); count])),
            &vec_scalar_mul(&n2, &z2)), &hadamard_multiply(&yn, &vec_scalar_mul(&sr, &u)));
        println!("Computed lu, ru");

        let [t0, t1, t2] = tx.access_tx();
        let t = t0 + (t1 * u) + (t2 * u * u);
        println!("Computed tu");

        let pi_lr = salt.alpha() + (salt.beta() * u);
        let pi_t = (z2 * salt.gamma()) + (salt.tau1() * u) + (salt.tau2() * u * u);
        println!("Computed all polynomial terms");

        let commit_c = inner_product(&l, &points.g_basis()) + inner_product(&r, &y_inv_h);
        assert_eq!(hadamard_multiply(&l, &r).iter().sum::<Scalar>(), t, "Prover system mess up!!");

        Self {
            l,
            r,
            t,
            pi_lr,
            pi_t,
            commit_c,
            y_inv_h,
        }
    }

    pub fn bullet_verifier(&self) -> (RistrettoPoint, [Scalar; 3]) {
        (self.commit_c, [self.t, self.pi_lr, self.pi_t])
    }

    pub fn lu(&self) -> Vec<Scalar> {
        self.l.clone()
    }

    pub fn ru(&self) -> Vec<Scalar> {
        self.r.clone()
    }

    pub fn tu(&self) -> Scalar {
        self.t
    }

    pub fn pi_lr(&self) -> Scalar {
        self.pi_lr
    }

    pub fn pi_t(&self) -> Scalar {
        self.pi_t
    }

    pub fn commit_c(&self) -> RistrettoPoint {
        self.commit_c
    }

    pub fn y_inv_h(&self) -> Vec<RistrettoPoint> {
        self.y_inv_h.clone()
    }
}

impl Salts {
    pub fn init() -> Self {
        let mut rng = OsRng;
        let alpha = Scalar::random(&mut rng);
        let beta = Scalar::random(&mut rng);
        let gamma = Scalar::random(&mut rng);
        let tau1 = Scalar::random(&mut rng);
        let tau2 = Scalar::random(&mut rng);

        Self {
            salt_alpha: alpha,
            salt_beta: beta,
            salt_gamma: gamma,
            salt_tau1: tau1,
            salt_tau2: tau2,
        }
    }

    pub fn alpha(&self) -> Scalar {
        self.salt_alpha
    }
    pub fn beta(&self) -> Scalar {
        self.salt_beta
    }
    pub fn gamma(&self) -> Scalar {
        self.salt_gamma
    }
    pub fn tau1(&self) -> Scalar {
        self.salt_tau1
    }
    pub fn tau2(&self) -> Scalar {
        self.salt_tau2
    }
}

impl ASVcommitment {
    pub fn compute(v: u64, range: usize, salt: &Salts, points: &GlobalPoints) -> Self {
        let salt_alpha = salt.alpha();
        let salt_beta = salt.beta();
        let salt_gamma = salt.gamma();

        let rng: u64 = rand::thread_rng().gen();
        let [mut al, mut ar] = ASVcommitment::compute_al_ar(v, range);
        let sl = gen_scalars(range, rng);
        let sr = gen_scalars(range, rng);

        println!("ASV commitment pre-requirements completed.");

        let commit_a_val = inner_product(&scalarize(&mut al), &mut points.g_basis())
            + inner_product(&scalarize(&mut ar), &mut points.h_basis()) + (salt_alpha * points.b_i());

        let commit_s_val = inner_product(&sl, &mut points.g_basis())
            + inner_product(&sr, &mut points.h_basis()) + (salt_beta * points.b_i());
        
        let commit_v_val = (Scalar::from(v) * points.g_i()) + (salt_gamma * points.b_i());

        println!("ASV commitments ready!");

        Self {
            commit_a: commit_a_val,
            commit_s: commit_s_val,
            commit_v: commit_v_val,
            al: scalarize(&mut al),
            ar: scalarize(&mut ar),
            sl,
            sr,
        }
    }

    fn compute_al_ar(v: u64, range: usize) -> [Vec<i64>; 2] {
        let mut al= to_bin(v.clone());
        while al.len() < range {
            al.push(0)
        }

        let ar = al.iter().map(|&x| x - 1).collect();

        [al, ar]
    }

    pub fn polynomial_const(&self) -> [Vec<Scalar>; 4] {
        [self.al.clone(), self.ar.clone(), self.sl.clone(), self.sr.clone()]
    }

    pub fn to_verifier(&self) -> [RistrettoPoint; 3] {
        println!("ASV commitments sent to verifier.");
        [self.commit_a, self.commit_s, self.commit_v]
    }
}