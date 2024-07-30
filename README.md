# ReVSeq Dashboard

This is a Django app for visualizing the ReVSeq data.

# Usage

To run a server:

```sh
python manage.py runserver
```

For Auspice integration:

```sh
nexstrain view dashboard/static -> changed to: npx auspice-revseq view --datasetDir dashboard/static/auspice
```

The home view will be at: [http://localhost:8000/home](http://localhost:8000/home).
