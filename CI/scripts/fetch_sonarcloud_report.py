import json
import os
import sys

import requests
import argparse

def parse_args():
    parser = argparse.ArgumentParser(description="Fetch SonarCloud metrics and save to a JSON file")
    parser.add_argument("token", type=str, help="SonarCloud token")
    parser.add_argument("project_key", type=str, help="SonarCloud project key")
    parser.add_argument("filename", type=str, help="Output filename")
    return parser.parse_args()

def export_sonar_metrics_to_json(token, project_key="acts-project_acts", filename="sonar_report.json"):
    # Endpoint for component measures (SonarCloud Web API)
    url = "https://sonarcloud.io/api/measures/component"

    metrics = "bugs,vulnerabilities,code_smells,coverage,duplicated_lines_density,security_rating,reliability_rating"

    params = {
        "component": project_key,
        "metricKeys": metrics,
        "branch": "main",
    }

    # 2. Fetch data from SonarCloud
    try:
        # Use Basic Auth: pass the token as the username with an empty password
        response = requests.get(url, params=params, auth=(token, ''))
        
        # Raise an exception for HTTP errors (4xx or 5xx)
        response.raise_for_status()
        
        # Parse the JSON response
        data = response.json()
        
        # 3. Save the result to a local JSON file
        # indent=4 makes the file human-readable (pretty-print)
        with open(filename, 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=4)
            
        print(f"Success! Report saved to: {filename}")

    except requests.exceptions.RequestException as e:
        print(f"An error occurred during the API request: {e}")
        raise SystemExit(1) from e
    except Exception as e:
        print(f"An unexpected error occurred: {e}")
        raise SystemExit(1) from e


if __name__ == "__main__":
 
    args = parse_args()
    token = args.token
    export_sonar_metrics_to_json(token, project_key="acts-project_acts", filename="sonar_report.json")