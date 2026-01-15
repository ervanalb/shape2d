# Shape2D Demo

Interactive egui-based demo for the Shape2D geometry library.

## Building for Desktop

```bash
cargo run
```

## Building for WebAssembly

To build for web deployment:

```bash
./build-wasm.sh
```

This will:
1. Compile the demo to WebAssembly
2. Generate JavaScript bindings using wasm-bindgen
3. Copy web assets to the `dist/` directory

## Testing Locally

After building, test the WASM version locally:

```bash
python3 -m http.server --directory dist 8080
```

Then open http://localhost:8080 in your browser.

## Deploying to GitHub Pages

The demo can be automatically deployed to GitHub Pages using the included GitHub Actions workflow. Simply push changes to the main branch, and the workflow will build and deploy the demo.

See `../.github/workflows/deploy.yml` for the workflow configuration.
