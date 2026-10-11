# Simulador cinemático: paquete desde un dron

Aplicación interactiva para el parcial de Física Mecánica.

## Desplegar en Streamlit Community Cloud

1. Crea un repositorio en GitHub.
2. Sube `app.py` y `requirements.txt` a la raíz del repositorio.
3. Entra a https://share.streamlit.io/ e inicia sesión con GitHub.
4. Selecciona **Create app**, elige tu repositorio, la rama y `app.py`.
5. Pulsa **Deploy**.

## Ejecutar localmente

```bash
pip install -r requirements.txt
streamlit run app.py
```

El modelo usa las ecuaciones por tramos del enunciado, con el eje vertical positivo hacia abajo.
