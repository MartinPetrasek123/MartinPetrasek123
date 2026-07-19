FROM python:3.12-slim

WORKDIR /analysis
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt
COPY . .
CMD ["python", "code/run_all.py"]
