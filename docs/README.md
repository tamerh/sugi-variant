# Sugi Variant

A **free, deterministic reference page for every clinically significant human variant** —
**463,882** pathogenic, likely-pathogenic and conflicting variants across **6,537** genes
(protein-coding, ncRNA and mitochondrial). Each page cross-sources ClinVar, gnomAD,
AlphaMissense, conservation, SpliceAI, ClinGen and more into one place, and — uniquely —
surfaces the **~109,000 variants where the computational evidence and the clinical call
disagree**.

Sugi Variant is the clinical-genetics reference member of the [sugi.bio](https://sugi.bio)
family — alongside **BioBTree** (identifier grounding), **Sugi Atlas** (gene/drug/disease
reference), **Sugi Predict** (patent-space target prediction) and **Enju** (workflow
orchestration). It resolves and grounds every identifier through BioBTree.

- **How it works** → [how-it-works.md](how-it-works.md) — the data, the resolver, and the
  interpretation rules.
- **Getting started** → [getting-started.md](getting-started.md) — searching, URLs, sets,
  and gene pages.

## What it is (and isn't)

- **Is:** a cross-sourced, always-cited *description* of the authoritative evidence for a
  variant — ClinVar's classification, ClinGen's applied ACMG codes, population frequency,
  in-silico predictors, mechanism — assembled deterministically so the same inputs always
  produce the same page.
- **Isn't:** a diagnostic device or a variant-classification engine. Sugi Variant **never
  computes its own classification, score, or verdict** — it reports what the sources say.
  Reference and research use only; not medical advice.

## The discordance layer

Every page reports whether the independent signals concur. The differentiator is the
inverse: the **evidence-disagreement** view collects the variants where they *don't* — a
predictor at odds with the clinical call, predictors split among themselves, or the
predictors unanimous where ClinVar's submitters conflict. It's a **quality-control signal
for review, never a reclassification**. See it per gene on any gene page.

## Free, with services

The public product is free and open (Google-crawlable). The business model is **services** —
custom builds on your own data and consultancy. For a custom build, <get in touch:
tamer.gur07@gmail.com>.
