# RC Column Seismic Model

Browser application for rectangular reinforced concrete columns. From the design data of a column it gives

- the drifts at yield, at peak force and at a 20 % loss of lateral strength (medians, 16–84 % bounds and design values exceeded by 95 % and 90 % of the tests), beside ASCE/SEI 41-17 and EN 1998-3;
- the shift of these medians when one or two tests of columns of the same type and programme are entered;
- the collapse fragility of a single-degree-of-freedom system;
- a Pinching4 material command for OpenSees (Tcl and OpenSeesPy).

The models were fitted to 299 cyclic tests from 64 test programmes, and their errors were obtained with each programme held out. The predictions apply within the range of these tests; outside it the application issues a warning.

**Use online:** https://USER.github.io/REPO/ (replace after publishing)

**Run locally:** open `index.html` in a browser. No installation or server is needed; all computations run in the browser.

## Files

| File | Content |
|---|---|
| `index.html` | user interface |
| `core.js` | computations (predictions, design values, test anchoring, code values, fragility, Pinching4 arguments) |
| `model.js` | fitted coefficients, reduction factors, cyclic-parameter medians and input ranges |

No test data are included.

## Citation

If you use this application, please cite the archived version (DOI on Zenodo) and the paper:

Tran, V.H., et al. Reliability of design-based predictions of nonlinear modelling parameters for reinforced concrete columns. Manuscript in preparation.

## License

MIT, see `LICENSE`.
