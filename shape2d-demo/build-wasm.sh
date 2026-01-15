#!/bin/bash
set -e

# Build the WASM binary
cargo build --release --target wasm32-unknown-unknown

# Install wasm-bindgen-cli if not already installed
# Note: Version must match the wasm-bindgen dependency in Cargo.toml
WASM_BINDGEN_VERSION="0.2.106"
if ! command -v wasm-bindgen &> /dev/null; then
    echo "Installing wasm-bindgen-cli v${WASM_BINDGEN_VERSION}..."
    cargo install wasm-bindgen-cli --version ${WASM_BINDGEN_VERSION}
else
    # Check if the installed version matches
    INSTALLED_VERSION=$(wasm-bindgen --version | cut -d' ' -f2)
    if [ "$INSTALLED_VERSION" != "$WASM_BINDGEN_VERSION" ]; then
        echo "Warning: Installed wasm-bindgen-cli version ($INSTALLED_VERSION) doesn't match required version ($WASM_BINDGEN_VERSION)"
        echo "Installing correct version..."
        cargo install -f wasm-bindgen-cli --version ${WASM_BINDGEN_VERSION}
    fi
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
