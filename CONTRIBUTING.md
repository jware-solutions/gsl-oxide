# Contributing

All contributions are welcome!

1. Fork the repository and create a branch for your changes.
2. Implement your changes and add or update tests when needed.
3. Run formatting, lint, and tests before opening a pull request:

   ```sh
   cargo fmt --all --check
   cargo clippy --workspace --all-targets --locked -- -D warnings
   cargo test --workspace --locked
   ```

4. Open a pull request describing your changes.

The repository and contributions are distributed under GPLv3 or later
(`GPL-3.0-or-later`); see [LICENSE](LICENSE).
