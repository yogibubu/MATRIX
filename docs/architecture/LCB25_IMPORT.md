# LCB25 Import

The LCB25 website provides a database of PCS2, SE and HPCS2 geometries. The
public page exposes three ZIP downloads:

- `PCS2.zip`
- `SE.zip`
- `HPCS2.zip`

ORACLE treats LCB25 as a remote geometry library source. The import flow is:

```text
LCB25 ZIP
  -> local library cache
  -> ORACLE-Babel import for each XYZ
  -> enriched XYZ with #SOURCE/#TOPOLOGY/#SYNTHONS
  -> fragment extraction or full-molecule reference search
```

LCB25 molecules can be used in two directions:

- as whole-molecule references for MORPHEUS/reference-assisted workflows;
- as fragment libraries after ORACLE topology/synthon preprocessing.

Conversely, an arbitrary query molecule can be fragmented by the same
`#TOPOLOGY/#SYNTHONS` state and compared against LCB25-derived fragments.

The adapter intentionally starts with URL planning and local archive extraction.
Search/index metadata should be added after the downloaded XYZ naming and any
sidecar metadata files are inspected.
