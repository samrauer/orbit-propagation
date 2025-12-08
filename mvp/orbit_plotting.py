import numpy as np
import plotly.graph_objects as go

from orbit_functions import R_EARTH, GaussSolution




def plot_gauss(sol: GaussSolution):

    fig = go.Figure()

    # earth Sphere
    theta = np.linspace(0, 2*np.pi, 50)
    phi = np.linspace(0, np.pi, 50)
    x_s = R_EARTH * np.outer(np.cos(theta), np.sin(phi))
    y_s = R_EARTH * np.outer(np.sin(theta), np.sin(phi))
    z_s = R_EARTH * np.outer(np.ones(50), np.cos(phi))
    fig.add_trace(go.Surface(x=x_s, y=y_s, z=z_s, colorscale="Blues", opacity=0.6, showscale=False, name="Earth"))


    # positions r1 and r2
    fig.add_trace(go.Scatter3d(
        x=[sol.r1[0]], y=[sol.r1[1]], z=[sol.r1[2]],
        mode="markers",
        marker=dict(
            size=6,
            color="red"
        ),
        name="r1"
    ))
    fig.add_trace(go.Scatter3d(
        x=[sol.r2[0]], y=[sol.r2[1]], z=[sol.r2[2]],
        mode="markers",
        marker=dict(
            size=6,
            color="green"
        ),
        name="r2"
    ))

    # velocity vectors (need to scale them to be visible)
    scale = 240

    # function to create a vector trace
    def plot_vector(r, v, color, name):
        end = r + v * scale
        return go.Scatter3d(
            x=[r[0], end[0]], y=[r[1], end[1]], z=[r[2], end[2]],
            mode="lines",
            line=dict(
                color=color, 
                width=6
            ),
            name=name
        )

    fig.add_trace(plot_vector(sol.r1, sol.v1, "red", "v1"))
    fig.add_trace(plot_vector(sol.r2, sol.v2, "green", "v2"))

    # layout settings
    fig.update_layout(
        title="gauss problem",
        scene=dict(
            xaxis_title="X (m)",
            yaxis_title="Y (m)",
            zaxis_title="Z (m)",
            aspectmode="data"
        ),
        autosize=False,
        width=600,
        height=400,
        template="plotly_dark"
    )

    return fig