KNU 500-trigger paired trial (2026-10-28)

Submit from the repository root:
  condor_submit condor/knu_500.sub
Do not run scripts/knu_500.sh directly on a login node.
Submitting again creates another trial using the same seed, not independent statistics.

One job reserves 1 CPU and 4 GiB RAM. It runs original, validates its closed
ROOT file, then runs light32 and validates it on the same worker node.
Each model collects 500 T1 AND T2 events (not 500 generated primaries).
Both use seed 261028500, distinct from the five-trigger seed 261026304.
Both reuse build-knu/job_2188629_0_bQsF6C and identical private patched
Geant4 11.4.1 Release libraries. Binary hashes and library paths are checked.
All baseline physics inputs remain unchanged; no new build is needed.

The wrapper and submit file intentionally pin this trial to target 500 and
seed 261028500. Change both consistently for a deliberately planned new trial.
Each simulation times out after 12 hours; the entire job is held after 25 hours.
A failed process or validation stops the pair and causes a Condor hold.
There is no automatic seed substitution. Inspect failures before retrying.

Outputs:
  runs/knu_500_20261028_1643/condor.CLUSTER.0.out
  runs/knu_500_20261028_1643/condor.CLUSTER.0.err
  runs/knu_500_20261028_1643/job_CLUSTER_0_RANDOM/
Each start gets a unique folder, preserving previous attempts.
Watch MODEL_START, MODEL_PASS, PAIR_PASS and WORKFLOW_EXIT_CODE in stdout.
The printed SIMULATION_LOG path contains detailed Geant4 output.
Final validation files: original_500/VALIDATION.json, light32_500/VALIDATION.json.
A successful pair writes summary/comparison_CLUSTER_0_job_*.json including
hit rates, binomial errors, path counts, optical totals, peak ADC values,
primary-event matching, paired errors where applicable, CPU/wall time and memory.
Primary equality does not imply optical-history equality.
One pair is insufficient to establish general speed or physics equivalence;
execution order, caching and changing host load remain possible confounders.

Monitor (replace CLUSTER with the submission number):
  condor_q CLUSTER -nobatch
  condor_q CLUSTER -af JobStatus RemoteHost
  tail -n 30 runs/knu_500_20261028_1643/condor.CLUSTER.0.out
  condor_history CLUSTER -limit 1 -af ClusterId JobStatus ExitCode
Completed status 4 and ExitCode 0, both validation PASS, and PAIR_PASS are
required before treating the pair as successful.

Preflight: shell syntax, Python compilation, 46-file bundle integrity,
prepared 500-trigger macros, Condor dry-run, and comparison regression on
the previously validated five-trigger pair all passed before submission.
