from fastapi import FastAPI
from pydantic import BaseModel
import pandas as pd
from shiny_app_prot import (
    expression_data_filtered,
    expression_data_clavata_filtered,
    predict_expression,
    generate_html_report,
)

app = FastAPI(title="Spider Gene API")

class SequenceRequest(BaseModel):
    sequence: str
    spider: str = "Bridge Spider"

@app.post("/predict")
async def predict(req: SequenceRequest):
    df = (
        expression_data_filtered
        if req.spider != "Clavata"
        else expression_data_clavata_filtered
    )
    tissues = df.columns[1:]
    preds = predict_expression(req.sequence, tissues)
    return {"tissues": list(tissues), "predicted_expression": preds.tolist()}

@app.get("/expression/{orthogroup}")
async def expression(orthogroup: str, spider: str = "Bridge Spider"):
    df = expression_data_filtered if spider != "Clavata" else expression_data_clavata_filtered
    row = df[df["Orthogroup"] == orthogroup]
    if row.empty:
        return {"error": "Orthogroup not found"}
    record = row.iloc[0].to_dict()
    return record

class ReportRequest(BaseModel):
    orthogroup: str
    spider: str = "Bridge Spider"
    sequence: str | None = None

@app.post("/report")
async def report(req: ReportRequest):
    df = expression_data_filtered if req.spider != "Clavata" else expression_data_clavata_filtered
    row = df[df["Orthogroup"] == req.orthogroup]
    if row.empty:
        return {"error": "Orthogroup not found"}
    expr_values = row.iloc[0, 1:].astype(float).values
    tissues = row.columns[1:]
    pred_vals = predict_expression(req.sequence, tissues) if req.sequence else None
    path = generate_html_report(req.orthogroup, tissues, expr_values, predicted_values=pred_vals)
    with open(path) as f:
        html = f.read()
    return {"report_html": html}
