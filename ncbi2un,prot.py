import requests
import time

def refseq_to_uniprot(refseq_id):
    # Step 1: Submit ID mapping job
    url = "https://rest.uniprot.org/idmapping/run"
    params = {
        "from": "RefSeq_Protein",
        "to": "UniProtKB",
        "ids": refseq_id
    }
    response = requests.post(url, data=params)
    response.raise_for_status()
    job_id = response.json()["jobId"]

    # Step 2: Poll for results
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        status_response = requests.get(status_url)
        status = status_response.json()
        if status.get("jobStatus") in ["RUNNING", "NEW"]:
            time.sleep(1)
        elif "results" in status or status.get("jobStatus") == "FINISHED":
            break
        else:
            raise RuntimeError("Mapping job failed or unexpected response")

    # Step 3: Get mapping results
    result_url = f"https://rest.uniprot.org/idmapping/uniprotkb/results/{job_id}"
    result_response = requests.get(result_url)
    result_response.raise_for_status()
    results = result_response.json()["results"]

    # Return all UniProt IDs
    return [entry["to"]["primaryAccession"] for entry in results]

# Example usage
refseq_id = "WP_001338221.1"
uniprot_ids = refseq_to_uniprot(refseq_id)
print(f"UniProt IDs for {refseq_id}: {', '.join(uniprot_ids)}")
