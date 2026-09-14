import logging
import hydra
from omegaconf import DictConfig
import torch
import scvi
import numpy as np
import random
from sklearn.model_selection import ParameterGrid
from utils.data_handling import load_dataset, find_datasets
import pipelines
import gc


@hydra.main(config_path="configs", config_name="config", version_base=None)
def main(cfg: DictConfig):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)

    ##### SETUP #####
    output_dir = hydra.core.hydra_config.HydraConfig.get().runtime.output_dir
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")

    # logging
    logger.info(f"Experimental results will be stored in: {output_dir}")
    logger.debug(f"Config: {cfg}")
    logger.info(f"Using device: {device}")

    # set random seed for reproducibility
    torch.manual_seed(cfg.seed)
    np.random.seed(cfg.seed)
    random.seed(cfg.seed)
    scvi.settings.seed = cfg.seed
    logger.info(f"Random seed set to: {cfg.seed}")

    ### datasets paths
    dat_paths = sorted(find_datasets(
        cfg.path,
        contains=cfg.contains,
        not_contains=["label"],
    ))
    assert len(dat_paths) > 0, "No datasets found"
    logger.info(f"Found {len(dat_paths)} datasets.")
    logger.debug(f"Datasets: {dat_paths}")

    ##### RUN EXPERIMENTS #####
    for i, p in enumerate(dat_paths):
        try:
            del dataset
            gc.collect()
            if device != "cpu":
                torch.cuda.empty_cache()
        except NameError:
            pass  # No dataset to delete in the first iteration

        dataset = load_dataset(p)

        logger.info(
            f"Run experiments on {len(cfg.pipeline)} pipelines: {list(cfg.pipeline.keys())}"
        )
        logger.info(f"Progress {(i+1)}/{len(dat_paths)}")
        logger.info(f"Dataset: {p}")
        logger.info(f"Dataset shape: {dataset.shape}")

        for pl_name, pl_params in cfg.pipeline.items():
            logger.info(f"Running pipeline: {pl_name}")

            param_grid = ParameterGrid(dict(pl_params))  # Create a grid of parameters
            logger.debug(f"Parameter grid: {param_grid}")
            logger.debug(f"Number of parameter combinations: {len(param_grid)}")
            for params in param_grid:
                logger.debug(f"Parameters: {params}")
                pl = getattr(pipelines, pl_name)(
                    **params, 
                    device=device, 
                    seed=cfg.seed,
                    unique_groups=cfg.unique_groups,
                )  # instantiate

                # Perform pipeline steps
                pl.preprocess(dataset)
                pl.fit()  # training
                # Log the device
                if isinstance(pl, pipelines.SCVI):
                    model_device = next(pl.model.module.parameters()).device
                    logger.info(f"Actual scVI model device: {model_device}")
                    if model_device.type == "cuda":
                        logger.info(f"GPU name: {torch.cuda.get_device_name(model_device)}")
                pl.get_latent()  # Embedding
                pl.visualization()  # TSNE/UMAP
                pl.clustering()  # Leiden
                pl.save(
                    f"{output_dir}/{pl_name}_{p.split('/')[-1].split('.')[0]}_{params}_{cfg.seed}",
                )

                # Clear GPU memory
                del pl
                gc.collect()
                if device != "cpu":
                    torch.cuda.empty_cache()

            logger.info(f"Pipeline {pl_name} completed for dataset {p}")


if __name__ == "__main__":
    main()
