# Contributing

All contributions are welcome!

1. Fork the repository and create a branch for your changes.
2. Implement your changes and add or update tests when needed.
3. Run the relevant checks listed below.
4. Open a pull request describing your changes.

The repository and contributions are distributed under GPLv3 or later
(`GPL-3.0-or-later`); see [LICENSE](LICENSE).

## Development commands

- Formatting: `cargo fmt --all --check`
- Lints: `cargo clippy --workspace --all-targets --locked -- -D warnings`
- Tests: `cargo test --workspace --locked`
- Documentation examples: `cargo test --workspace --doc --locked`
- Benchmarks: `cargo bench --workspace --locked`
- API documentation: `cargo doc --workspace --no-deps --locked`
- Build one crate: `cargo build -p <package> --locked`

Use `-p <package>` instead of `--workspace` to check a single crate.

## References and testing

The public [GSL 2.8 definitions and documented behavior](https://www.gnu.org/software/gsl/doc/html/) are the behavior reference. Short mathematical provenance notes live beside the Rust source.
