#!/bin/bash
set -e

# Build the WASM binary
cargo build --release --target wasm32-unknown-unknown

# Install wasm-bindgen-cli if not already installed
if ! command -v wasm-bindgen &> /dev/null; then
    echo "Installing wasm-bindgen-cli..."
    cargo install wasm-bindgen-cli
fi

# Generate JavaScript bindings
wasm-bindgen --out-dir ./dist \
    --target web \
    --no-typescript \
    ../target/wasm32-unknown-unknown/release/shape2d-demo.wasm

# Copy assets
cp assets/index.html dist/

echo "Build complete! Output is in shape2d-demo/dist/"
echo "To test locally, run: python3 -m http.server --directory dist 8080"
echo "Then open http://localhost:8080 in your browser"
