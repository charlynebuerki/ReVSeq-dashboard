# Use the official Python image from the Docker hub
FROM python:3.11-slim

# Set the working directory in the container
WORKDIR /app

# Copy the requirements.txt file to install the dependencies
COPY requirements.txt .

RUN pip install --no-cache-dir -r requirements.txt

# Copy the entire Django project to the container
COPY . .

# Expose the port the Django app runs on
EXPOSE 8000

# Set environment variables (optional)
#ENV DJANGO_SETTINGS_MODULE=myproject.settings
#ENV PYTHONUNBUFFERED=1

# Run the Django dashboard using your custom command
CMD ["python", "manage.py", "runserver", "0.0.0.0:8000" ]
