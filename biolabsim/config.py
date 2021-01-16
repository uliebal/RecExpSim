
from pathlib import Path


# Directories will be resolved relative to the position of this file.
_rootpath = Path(os.path.abspath(__file__)).parent.parent


# Directory containing all biological models.
MODEL_DIR: Path = _rootpath / "data" / "metabolic_model"
