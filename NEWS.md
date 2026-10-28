# NECorr 3.6.0

* The default ranking is now the rank product of up-regulation and the number
  of targets of a gene (`degree_type = "targets"`) for directed networks.
  Undirected networks keep the number of partners. One rule now serves the
  human benchmark and the plant applications: on human CRISPR dependencies it
  equals total degree among the first ten genes on both networks and is 0.7
  hits behind among the first 50 on one of them; in rice it places about twice
  as many characterised genes among the first 50.
* Counting down-regulation together with up-regulation (`direction = "both"`)
  was tested with the number of targets, the number of regulators and total
  degree. It lowers the number of hits on the human benchmark in every form
  and stays an option, not the default.

# NECorr 3.5.0

* New argument `degree_type`. `"total"` (default, unchanged) counts all partners
  of a gene. `"targets"` counts only its targets (out-degree) in a directed
  network. On human CRISPR dependencies most of the gain of degree comes from
  the number of targets: the two settings give the same number of hits among
  the first ten genes on both networks and at every list length on CollecTRI,
  and total degree keeps an advantage among the first 50 and 100 genes on
  DoRothEA. With a rice network predicted from binding motifs (PlantRegMap),
  the number of regulators predicted for a gene lowers the ranking, and
  `"targets"` places about twice as many characterised genes among the first 50
  as up-regulation alone. That result rests on one crop.

# NECorr 3.4.0

* New default ranking: rank product of up-regulation in the target condition (`DE` with `direction = "up"`) and `Degree`
  (`criteria = c("DE", "Degree")`, `scaling = "rank"`, `aggregation = "product"`, `weights = "equal"`).
  On DepMap 26Q1 it recovers more lineage-selective CRISPR dependencies among the first 50 genes than the six-criteria
  scheme of 3.3, than network propagation and than expression alone, and a supervised model given all criteria does not improve on it.
* New arguments `direction` ("up", "both") and `scaling` ("rank", "softcap").
* `rank_regulators()` orders candidates by their own score (`order_by = "score"`); the target-enrichment order of earlier versions is `order_by = "enrichment"`.
* The 3.3 scheme remains available: `criteria = c("Coexpr", "CondCoexpr", "TS", "DE", "Bridging", "Degree")`, `weights = "ahp"`,
  `aggregation = "sum"`, `direction = "both"`, `scaling = "softcap"`.

