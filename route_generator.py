from planner_engine import plan_multistep


def generate_routes(smiles: str, max_steps: int, n_routes: int = 3):
    routes = []

    for i in range(n_routes):
        plan = plan_multistep(smiles, max_steps)

        route = {
            "route_id": i + 1,
            "strategy": f"strategy_{i+1}",
            "steps": plan["steps"],
            "overall_confidence": plan["overall_confidence"]
        }

        routes.append(route)

    return routes
