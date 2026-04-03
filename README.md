## Usage

### Option 1: Command Line

1. Prepare your `sample_table.csv` in the data directory.
2. Run the pipeline:
```bash
sbatch -o slurms_out.log nanopore_pipeline.sh -b 20260209 -i /path/to/data -o /path/to/results
```

### Option 2: Web UI

1. Launch the Streamlit app:
```bash
streamlit run src/app.py
```

2. Fill in the parameters in the sidebar:
   - **Batch ID**: your batch name (e.g. `20260209`)
   - **Data Dir**: path to your raw data directory (must contain `sample_table.csv`)
   - **Results Dir**: path to store the output
   - **Species**: `Mouse` or `Human`

3. Select run mode:
   - **Local**: runs directly on the current machine
   - **Slurm**: submits to the cluster via `sbatch`, and monitors progress in real time

4. Click **🚀 Run MiXCR pipeline** to start.

5. After the pipeline finishes, go to **📊 Results** to explore:
   - **🔍 Read Analysis**: read retention summary and per-step filtering stats
   - **📊 Clone Analysis**: clone table, SHM distribution, CDR3 length distribution, V gene usage