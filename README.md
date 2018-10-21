# actin-entropy-paper
Repository for code used in "Entropy Production Rate is Maximized in Non-contractile Actomyosin"

### Entropy production
Depending on type of data to analyze (JFilament tracking data, cytosim outputs, or axoneme data), use the corresponding `save____Output.m` file to extract mode amplitude time series and filament lengths. Feed the resulting time series into `stochasticEntropyChange.m` to get entropy change time series.
