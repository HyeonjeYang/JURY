# Example workflow

Use your own map, insert, and reads stored outside this repository:

```sh
jury verify /path/to/assembled.dna /path/to/reads --insert /path/to/gblock.dna --orf-motif ATGGGATCT
jury primers /path/to/assembled.dna --spacing 600
```

Read the insert/junction results, feature verdicts, and coverage together.
`LOW-CONFIDENCE` means the jury is still out, even when several reads agree.
A `MISMATCH` in a low-quality stretch needs trace review. Primer walk coverage
is a prediction; new primers still need experimental validation.

No experimental sequences are distributed with JURY.
