### v0.4.2 — 2026-04-11
- **Fixed:** bug of `+I` model fixed (in v0.4.0 and v0.4.1 it could lead to NA for lnL thus no output for the Hessian).

### v0.4.1 — 2026-03-27
- **Fixed:** bug of not uploading some dependencies in the folder `additional_scripts/`

### v0.4.0 — 2026-03-27
- **Fixed:** bug fixed for LG4M combined with `+I`.
- **Improved:** speed up for lnL calculation via optimization of the matrix exponential.
- **Improved:** speed up for Hessian calculation particularly under the 2nd-order derivative of lnL (`--hessian_type fd`).

