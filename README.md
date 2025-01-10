
# Zero-Knowledge Range Proof with Bulletproof + Linear Verification

Implementation of **Zero-Knowledge Range Proof** for ranges of the form `2^n`. The system uses **Bulletproofs O(log n)** for the default verification method, but **Linear Verification O(n)** can also be used (currently commented out). The range proofs ensure that the prover knows a value within a specified range without revealing any other information about the value itself.

## Features

- **Zero-Knowledge Proofs**: The proof ensures that the prover knows a value within the range `2^n` without revealing the actual value, preserving privacy.
- **Bulletproof Verification**: Bulletproofs are used for efficient and compact range proofs, where the verifier can be confident that the prover knows a valid value within the range.
- **Linear Verification**: An optional linear verification method is provided (commented out), offering an alternative verification approach.
- **Custom Prover and Verifier System**: A custom prover and verifier system is designed with specifically allotted structs and implementations for the range proof.
- **Efficient Implementation**: The main file contains a step-by-step process for generating commitments, computing proofs, and verifying them.

## Mechanism

### 1. **Prover Setup**:
   - The prover selects a value `v` and converts it to a binary vector `al`. 
   - The prover then computes commitments using the `ASVcommitment` and `T1T2commitment` structs, and performs the necessary computations for Bulletproof verification.

### 2. **Global Points Generation**:
   - Global points are generated using a specific `seed` and are used to define the basis vectors `g_basis` and `h_basis`. 
   - A random point `g_i` and `b_i` are also generated for use in the range proof.

### 3. **Commitment Generation**:
   - The prover generates a series of commitments for the value `v`, including the `T1` and `T2` commitments, which are used for Bulletproof verification.

### 4. **Bulletproof Verification** :
   - Bulletproofs are used to verify that the prover knows a value within the range `2^n`. This is achieved through a series of scalar and point operations, which are checked by the verifier to ensure that the proof is valid.
   - Time complexity: O(log n)

### 5. **Linear Verification**:
   - (Currently commented out) An optional linear verification method is also available, which performs a more direct check of the range proof.
   - Time complexity: O(n)

## Dependencies

The following dependencies are required to build and run this project:

- `curve25519-dalek`: A Rust library implementing the Curve25519 elliptic curve and the Ristretto group.
- `rand`: A random number generation library used for cryptographic randomness.
- `rand_core`: Core random number generation functionality for Rust.
- `rand_chacha`: A Chacha random number generator, used for deterministic random number generation with a seed.
- `sha2`: A Rust implementation of the SHA-2 cryptographic hash function.

### Add the following to your `Cargo.toml`:

```toml
[dependencies]
curve25519-dalek = {version = "4.*", features = ["rand_core"]}
rand = "0.8"
rand_core = "*"
rand_chacha = "0.3"
sha2 = "0.10"
```

## Steps to Run

1. **Clone the Repository**:
   ```bash
   git clone https://github.com/nikillxh/zk-range-proof
   cd zk-range-proof
   ```

2. **Build and Run**:
   To compile and run the code:
   ```bash
   cargo build
   cargo run
   ```

3. **Configuration**:
   - The range can be modified by changing the `range` value (which represents `2^n`).
   - You can also toggle between Bulletproof and Linear verification by uncommenting or commenting the relevant lines in the `main.rs` file.

## Files and Modules

- `main.rs`: The entry point of the program that coordinates the prover and verifier steps.
- `prover.rs`: Contains logic for computing commitments and performing the proof.
- `verifier.rs`: Handles the verification process for Bulletproof and linear verification.
- `generator.rs`: Defines global points generation and scalar creation.
- `operations.rs`: Contains utility functions for operations like vector addition, scalar folding, etc.
- `bulletproof.rs`: Implements the Bulletproof protocol for range proof verification.

## References
- **Rareskills ZK-Book**: https://www.rareskills.io/zk-book
- **Bulletproofs**: https://eprint.iacr.org/2017/1066.pdf

## License

This project is licensed under the GPL-3.0 License - see the [LICENSE](LICENSE) file for details.
