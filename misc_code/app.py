import kaleido
import plotly.express as px
import plotly.io as pio

# Create a simple plot
fig = px.scatter(x=[1, 2, 3], y=[4, 5, 6], title="Test Plot")

# Try saving as a PNG
try:
    pio.write_image(fig, "kaleido_test.png")
    print("Kaleido test successful. Plot saved as 'kaleido_test.png'.")
except Exception as e:
    print(f"Kaleido test failed: {e}")
