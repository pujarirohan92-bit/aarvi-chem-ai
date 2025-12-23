from fastapi import FastAPI
from pydantic import BaseModel
from chemistry_router import route_by_chemistry
from route_generator import generate_routes
from route_ranker import rank_routes
from inorganic_engine import plan_inorganic

app = FastAPI(title="Aarvi Chem AI – V2.2 (Organic + Inorganic)")


class PredictRequest(BaseModel):
    smiles: str
    max_steps: int = 3
    n_routes: int = 3


@app.post("/predict")
def predict(req: PredictRequest):
    routing = route_by_chemistry(req.smiles)
    engine = routing["engine"]

    # ORGANIC (small + large)
    if engine in ["organic_engine", "organic_large_engine"]:
        routes = generate_routes(req.smiles, req.max_steps, req.n_routes)
        ranked = rank_routes(routes)

        return {
            "target_smiles": req.smiles,
            "version": "v2.2",
            "summary": {
                "best_route_id": ranked[0]["route_id"],
                "overall_confidence": ranked[0]["overall_confidence"],
                "overall_expected_yield": ranked[0].get("overall_expected_yield")
            },
            "routes": ranked,
            "meta": {
                "chemistry_type": routing["info"]["type"],
                "note": routing["info"]["note"]
            }
        }

    # INORGANIC / COORDINATION
    if engine == "inorganic_engine":
        route = plan_inorganic(req.smiles, req.max_steps)

        return {
            "target_smiles": req.smiles,
            "version": "v2.2",
            "summary": {
                "best_route_id": route["route_id"],
                "overall_confidence": route["overall_confidence"]
            },
            "routes": [route],
            "meta": {
                "chemistry_type": routing["info"]["type"],
                "note": "Coordination chemistry baseline (explainable)"
            }
        }

    # ORGANOMETALLIC / UNKNOWN
    return {
        "target_smiles": req.smiles,
        "version": "v2.2",
        "summary": None,
        "routes": [],
        "meta": {
            "chemistry_type": routing["info"]["type"],
            "note": "Detected but engine under development"
        }
    }
