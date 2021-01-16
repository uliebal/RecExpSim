
from pathlib import Path


# Directories will be resolved relative to the position of this file.
_rootpath = Path(os.path.abspath(__file__)).parent.parent


# Directory containing all biological models.
METABOLIC_MODEL_DIR: Path = _rootpath / "data" / "metabolic_model"


# Directory containing all expression predictors.
EXPRESSION_PREDICTOR_DIR: Path = _rootpath / "data" / "expression_predictor"
