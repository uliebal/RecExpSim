
from FermProSimFun import Monod_Model


for i in range(2):
	TestModel = Monod_Model()
	model_results = TestModel.calculate_monod()
	TestModel.plot_results()
	TestModel.to_json()
