# Tool page header styling
header_css = """
<style>
:root {
    --color: #333;
}

div.header-container {
    position: relative;
    left: 50%;
    right: 50%;
    margin-left: -50vw;
    margin-right: -50vw;
    margin-top: -9rem;
    width: 100vw;
    min-height: 300px;
    box-sizing: border-box;
    text-align: center;
    color: var(--color);
    display: flex;
    flex-direction: column;
    align-items: center;
    justify-content: center;
    background-size: cover;
    background-position: center 90%;
    background-attachment: fixed;
    overflow: hidden;
    margin-bottom: 2rem;
}

div.header-container::before {
    content: "";
    position: absolute;
    top: 0;
    left: 0;
    right: 0;
    bottom: 0;
    backdrop-filter: blur(0px); /* gaussian blur */
    z-index: 1;
    pointer-events: none;
}

div.header-container img {
    z-index: 3;
    position: relative;
}

@media screen and (max-width: 799px) {
    div.header-container {
        min-height: 250px;
    }
}
</style>
"""

bubble_css = """
<style>
.success-bubble,
.info-bubble,
.warning-bubble,
.error-bubble {
    display: flex;
    flex-direction: column;
    align-items: flex-start;
    padding: 1rem;
    border-radius: 0.5rem;
    margin: 1rem 0;
    transition: all 0.3s ease-in-out;
    line-height: 1.5;
    font-size: 1.125rem;
    font-weight: 300;
    white-space: normal;
}

.success-bubble b,
.info-bubble b,
.warning-bubble b,
.error-bubble b {
    display: inline;
    font-weight: bold;
    margin: 0 !important;
    padding: 0 !important;
    white-space: normal;
}

.success-bubble p,
.info-bubble p,
.warning-bubble p,
.error-bubble p {
    display: inline;
    margin: 0 !important;
}

.success-bubble {
    background-color: #E8F5E9;
    border-left: 4px solid #22C55E;
    color: #14532D;
}
.success-bubble:hover {
    background-color: #D1FAE5;
    transform: scale(1.05);
}

.info-bubble {
    background-color: #DBEAFE;
    border-left: 4px solid #3B82F6;
    color: #1E3A8A;
}
.info-bubble:hover {
    background-color: #BFDBFE;
    transform: scale(1.05);
}

.warning-bubble {
    background-color: #FEF9C3;
    border-left: 4px solid #EAB308;
    color: #713F12;
}
.warning-bubble:hover {
    background-color: #FEF08A;
    transform: scale(1.05);
}

.error-bubble {
    background-color: #FEE2E2;
    border-left: 4px solid #EF4444;
    color: #7F1D1D;
}
.error-bubble:hover {
    background-color: #FECACA;
    transform: scale(1.05);
}
</style>
"""

button_css = """
<style>
.stButton>button[kind="primary"],
.stFormSubmitButton>button,
[data-testid="stSidebarContent"] [data-testid="stPopoverButton"] {
    text-transform: uppercase;
    text-decoration: none;
    padding: 7px 40px;
    display: inline-flex;
    align-items: center;
    border-radius: 100px;
    transition: all .2s;
    margin: 0 auto;
    font-size: 20px;
    position: relative;
    background-color: #dc474b;
    color: #fff;
    border: 2px solid #dc474b;
}

.stButton>button[kind="primary"]>div>svg,
.stFormSubmitButton>button>div>svg {
    margin-right: 8px;
}

.stButton>button[kind="primary"]:hover,
.stFormSubmitButton>button:hover,
[data-testid="stSidebarContent"] [data-testid="stPopoverButton"]:hover {
    transform: translateY(-0.5px);
    box-shadow: 0 5px 10px rgba(0, 0, 0, 0.1);
    background-color: #e55a5e;
    color: #fff;
}

.stButton>button[kind="primary"]::before,
.stFormSubmitButton>button::before,
[data-testid="stSidebarContent"] [data-testid="stPopoverButton"]::before {
    content: "";
    width: 100%;
    height: 100%;
    inset: 0;
    position: absolute;
    background-color: transparent;
    border-radius: 100px;
    transition: all 0.3s ease;
}

.stButton>button[kind="primary"]:hover::before,
.stFormSubmitButton>button:hover::before,
[data-testid="stSidebarContent"] [data-testid="stPopoverButton"]:hover::before {
    box-shadow: 0 0 15px 5px #dc474b;
    color: #fff;
}

.stButton>button[kind="primary"]:active,
.stFormSubmitButton>button:active,
[data-testid="stSidebarContent"] [data-testid="stPopoverButton"]:active {
    transform: translateY(0px);
    box-shadow: 0 3px 5px rgba(0, 0, 0, 0.1);
}

.stButton>button[kind="primary"]::after,
.stFormSubmitButton>button::after,
[data-testid="stSidebarContent"] [data-testid="stPopoverButton"]::after {
    content: "";
    display: inline-block;
    height: 100%;
    width: 100%;
    border-radius: 100px;
    position: absolute;
    top: 0;
    left: 0;
    z-index: -1;
    transition: all .4s;
    background-color: #dc474b;
}

.stButton>button[kind="primary"]:hover::after,
.stFormSubmitButton>button:hover::after,
[data-testid="stSidebarContent"] [data-testid="stPopoverButton"]:hover::after {
    transform: scaleX(1.2) scaleY(1.3);
    opacity: 0;
}

.stButton>button[kind="secondary"] {
    text-transform: uppercase;
    text-decoration: none;
    padding: 7px 40px;
    display: inline-flex;
    align-items: center;
    border-radius: 100px;
    margin: 0 auto;
    position: relative;
    background-color: #fff;
    color: #000;
    border: 2px solid #dc474b;
    transition: all 0.3s ease;
}

.stButton>button[kind="secondary"]>div>svg {
    margin-right: 8px;
}

.stButton>button[kind="secondary"]:hover {
    background-color: #dc474b;
    color: #fff;
    transition: all 0.3s ease;
}

.stButton>button[kind="secondary"]::before {
    content: "";
    width: 100%;
    height: 100%;
    inset: 0;
    position: absolute;
    background-color: transparent;
    border-radius: 100px;
    transition: all 0.3s ease;
}

.stButton>button[kind="secondary"]:hover::before {
    box-shadow: 0 0 10px 3px #dc474b;
}

[data-testid="stFormSubmitButton"] {
    color: #fff !important;
}

.stButton>button[kind="tertiary"] {
    text-decoration: none;
    padding: 7px 20px;
    display: inline-flex;
    align-items: center;
    margin: 0 auto;
    font-size: 18px;
    background-color: #f0f2f6;
    color: #000;
    border: 1.5px solid #d5d6d8;
    border-radius: 20px;
    box-shadow: 0 1px 2px rgba(0,0,0,0.10);
    transition: background 0.2s, color 0.2s, border 0.2s;
    cursor: pointer;
}

.stButton>button[kind="tertiary"]>div>svg {
    margin-right: 8px;
}

.stButton>button[kind="tertiary"]:hover {
    background-color: #fff;
    color: #dc474b;
    border: 1.5px solid #dc474b;
    transition: all 0.3s ease;
}

div[data-testid="stColumn"]:has(.stButton) {
    width: fit-content !important;
    flex: unset;
}

div[data-testid="stColumn"]:has(.stButton) .stButton {
    width: fit-content !important;
}

/* Force text in form submit button to remain white after clicked */
[data-testid="stFormSubmitButton"] [data-testid="stBaseButton-secondaryFormSubmit"]:focus:not(:active),
[data-testid="stFormSubmitButton"] [data-testid="stBaseButton-secondaryFormSubmit"]:focus:not(:active) *,
[data-testid="stFormSubmitButton"] [data-testid="stBaseButton-secondaryFormSubmit"]:focus:not(:active) p,
[data-testid="stFormSubmitButton"] [data-testid="stBaseButton-secondaryFormSubmit"]:focus:not(:active) div,
[data-testid="stFormSubmitButton"] [data-testid="stBaseButton-secondaryFormSubmit"]:focus:not(:active) span {
    color: #fff !important;
}
</style>
"""

modify_streamlit_style = """
<style>
/* ──────────────────────────────
        Toolbar (Top)
──────────────────────────────── */

/* Hide streamlit toolbar above the page */
div[data-testid="stToolbar"] {
    visibility: hidden;
    height: 0%;
    position: fixed;
}

div[data-testid="stDecoration"] {
    visibility: hidden;
    height: 0%;
    position: fixed;
}

div[data-testid="stStatusWidget"] {
    visibility: hidden;
    height: 0%;
    position: fixed;
}

#MainMenu {
    visibility: hidden;
    height: 0%;
}

header {
    visibility: hidden;
    height: 0%;
}

footer {
    visibility: hidden;
    height: 0%;
}

/* Control the padding of the main content, usually padding top should be 3rem */
.block-container {
    padding-top: 0rem !important;
    padding-bottom: 5rem;
    padding-left: 5rem;
    padding-right: 5rem;
}

[data-testid="stBaseButton-headerNoPadding"] * {
    width: 1.2em;
    height: 1.2em;
    color: #dc474b !important;
}

/* Allow sidebar toggle button to remain visible */
[data-testid="stHeaderActionElements"] button:not([kind="header"]) {
    display: none !important;
}

/* Ensure sidebar collapse button is visible */
button[kind="header"] {
    display: block !important;
    visibility: visible !important;
}

/* ──────────────────────────────
        Sidebar (Left)
──────────────────────────────── */

[data-testid="stSidebarContent"] {
    background: #0e1117 !important;
    color: #fff !important;
}

[data-testid="stSidebarContent"] * {
    color: #fff !important;
}

[data-testid="stSidebarContent"] [data-testid="stBaseButton-secondary"] *,
[data-testid="stSidebarContent"] [data-testid="stTextInputRootElement"] *,
[data-testid="stSidebarContent"] [data-testid="stSelectbox"] * {
    color: initial !important;
}

[data-testid="stSidebarContent"] [data-testid="stSelectbox"] [data-testid="stMarkdownContainer"] * {
    color: #fff !important;
}

[data-testid="stSidebarContent"] [data-testid="stExpander"] {
    border: 1.5px solid #3F3F3F !important;
    border-radius: 8px !important;
    overflow: hidden;
    background: #232323 !important;
}

[data-testid="stSidebarContent"] [data-testid="stExpander"] a:visited {
    background: #fff !important;
}

[data-testid="stSidebarContent"] .stButton>button[kind="tertiary"] {
    color: #fff !important;
    background-color: #232323 !important;
    border: 1.5px solid #3F3F3F !important;
}

[data-testid="stSidebarContent"] .stButton>button[kind="tertiary"]:hover {
    background-color: #3F3F3F !important;
    border: 1.5px solid #dc474b !important;
}

[data-testid="stSidebarContent"] .info-bubble {
    color: #000 !important;
}

div[data-baseweb="popover"][data-testid="stPopoverBody"] {
    min-width: 600px !important;
    max-width: 700px !important;
}

[data-testid="stBaseButton-elementToolbar"], 
[data-testid="stElementToolbarButtonIcon"] {
    fill: #fff !important;
    background-color: #ec3d41 !important;
}

[data-testid="stSpinner"], 
[data-testid="stSpinner"] > div {
    display: flex;
    justify-content: center;
    align-items: center;
}

[data-testid="stToastContainer"] {
    font-size: 10rem;
}

[data-testid="stHeader"] {
    height: 0rem;
}

[data-testid="stSidebarContent"] [data-testid="stNumberInputContainer"] * {
    color: #000 !important;
}

[data-testid="stSidebarContent"] [data-testid="stFileUploaderDropzone"] * {
    color: #000 !important;
}
</style>
"""

custom_metric_delta_css = """
<style>
[data-testid="stMetricDelta"] svg {
    display: none !important;
}
</style>
"""

dialog_and_plot_css = """
<style>
div[data-testid="stDialog"] div[role="dialog"] {
    width: 70% !important;
    max-width: none;
    background: rgba(255, 255, 255, 0.85);
    border-radius: 16px;
    border: 1px solid rgba(255, 255, 255, 0.4);
    backdrop-filter: blur(10px);
    -webkit-backdrop-filter: blur(10px);
    box-shadow: 0 8px 32px rgba(0, 0, 0, 0.15);
    color: #000000 !important;
}

div[data-testid="stDialog"] div[role="dialog"] * {
    color: #000000 !important;
}

div[data-testid="stDialog"] div[role="dialog"] h1,
div[data-testid="stDialog"] div[role="dialog"] h2,
div[data-testid="stDialog"] div[role="dialog"] h3,
div[data-testid="stDialog"] div[role="dialog"] h4 {
    margin-top: 0 !important;
    padding-top: 0 !important;
    color: #000000 !important;
}

div[data-testid="stDialog"] div[role="dialog"] a {
    color: #000000 !important; 
    text-decoration: underline;
}

div[data-testid="stDialog"] div[role="dialog"] a:hover {
    color: #333333 !important;
}

div[data-testid="stDialog"] div[role="dialog"] div[data-testid="stHorizontalBlock"] {
    display: flex !important;
    align-items: center !important;
}

div[data-testid="stDialog"] div[role="dialog"] div[data-testid="stHorizontalBlock"] > div {
    display: flex !important;
    align-items: center !important;
    min-height: 38px !important;
}

div[data-testid="stDialog"] div[role="dialog"] div[data-testid="stButton"] {
    display: flex !important;
    justify-content: center !important;
    align-items: center !important;
    width: 100% !important;
}

div[data-testid="stDialog"] div[role="dialog"] div[data-testid="stButton"] button {
    display: flex !important;
    justify-content: center !important;
    align-items: center !important;
    width: 100% !important;
    height: 100% !important;
    padding: 0.5rem 1rem !important;
}

div[data-testid="stDialog"] div[role="dialog"] div[data-testid="stButton"] button p {
    margin: 0 !important;
    padding: 0 !important;
    display: flex !important;
    align-items: center !important;
    justify-content: center !important;
    gap: 0.5rem !important;
}

.js-plotly-plot .plotly .main-svg text,
.js-plotly-plot .plotly .main-svg .gtitle,
.js-plotly-plot .plotly .main-svg .xtick text,
.js-plotly-plot .plotly .main-svg .ytick text,
.js-plotly-plot .plotly .main-svg .legend text,
.js-plotly-plot .plotly .main-svg .annotation-text,
.js-plotly-plot .plotly .main-svg .infolayer text {
    fill: #000000 !important;
    color: #000000 !important;
}

.js-plotly-plot .plotly .main-svg .gtitle text,
.js-plotly-plot .plotly .main-svg .g-gtitle text,
.js-plotly-plot .plotly .main-svg .xaxis-title text,
.js-plotly-plot .plotly .main-svg .yaxis-title text {
    fill: #000000 !important;
    color: #000000 !important;
}

.js-plotly-plot .plotly .main-svg .hovertext {
    color: #000000 !important;
    background-color: rgba(255, 255, 255, 0.9) !important;
}
</style>
"""

metrics_tiles_hover_css = """
<style>
.metrics-tile-container {
    position: relative;
    width: 100%;
    height: auto;
    min-height: 160px;
    margin: 8px 0;
    perspective: 1000px;
    cursor: pointer;
    animation: tile-entrance 0.5s ease-out;
}

.metrics-tile-front {
    position: relative;
    width: 100%;
    min-height: 160px;
    backface-visibility: hidden;
    border-radius: 12px;
    padding: 10px;
    color: #333;
    text-align: center;
    box-shadow: 0 6px 20px rgba(0,0,0,0.1);
    border: 1px solid #d3d3d3;
    display: flex;
    flex-direction: column;
    justify-content: center;
    transition: transform 0.6s ease-in-out;
    transform-style: preserve-3d;
    background: inherit;
    z-index: 2;
    overflow: hidden;
}

/* Modern gradient backgrounds with circular elements */
.metrics-tile-front::before {
    content: "";
    position: absolute;
    top: -50px;
    right: -50px;
    width: 100px;
    height: 100px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-front::after {
    content: "";
    position: absolute;
    bottom: -30px;
    left: -30px;
    width: 60px;
    height: 60px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-front h3 {
    line-height: 1.0 !important;
    padding: 0 !important;
}

.metrics-tile-front h3 + div {
    line-height: 1.0 !important;
    padding: 0 !important;
}

.metrics-tile-back {
    position: absolute;
    top: 0;
    left: 0;
    width: 100%;
    height: 100%;
    backface-visibility: hidden;
    border-radius: 12px;
    padding: 15px;
    color: #000;
    text-align: center;
    display: flex;
    flex-direction: column;
    justify-content: center;
    align-items: center;
    transform: rotateX(180deg);
    transition: transform 0.6s ease-in-out;
    transform-style: preserve-3d;
    font-size: 16px;
    line-height: 1.5;
    word-wrap: break-word;
    overflow: hidden;
    font-weight: 500;
    box-shadow: 0 6px 20px rgba(0,0,0,0.1);
    border: 1px solid #d3d3d3;
    z-index: 1;
}

.metrics-tile-back h4 {
    margin: 0 0 8px 0;
    font-size: 21px;
    font-weight: 600;
    line-height: 1.0 !important;
    color: inherit;
}

.metrics-tile-back div {
    font-size: 17px;
    line-height: 1.0 !important;
    color: inherit;
}

.metrics-tile-container:hover .metrics-tile-front {
    transform: rotateX(180deg);
    animation: pulse-glow 2s infinite;
}

.metrics-tile-container:hover .metrics-tile-back {
    transform: rotateX(360deg);
}

/* Modern dark gradients for different metric types */
.metrics-tile-back.ic50-bg {
    background: linear-gradient(135deg, #1a1a2e 0%, #16213e 50%, #0f3460 100%);
    color: #fff;
}

.metrics-tile-back.ic50-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.ic50-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.pic50-bg {
    background: linear-gradient(135deg, #2d1b69 0%, #11998e 50%, #38ef7d 100%);
    color: #fff;
}

.metrics-tile-back.pic50-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.pic50-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.affinity-bg {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: #fff;
}

.metrics-tile-back.affinity-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.affinity-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.confidence-bg {
    background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
    color: #fff;
}

.metrics-tile-back.confidence-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.confidence-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.similarity-bg {
    background: linear-gradient(135deg, #fa709a 0%, #fee140 100%);
    color: #333;
}

.metrics-tile-back.similarity-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.15);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.similarity-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.12);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.plddt-bg {
    background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%);
    color: #333;
}

.metrics-tile-back.plddt-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.15);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.plddt-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.12);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.ptm-bg {
    background: linear-gradient(135deg, #ffecd2 0%, #fcb69f 100%);
    color: #333;
}

.metrics-tile-back.ptm-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.15);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.ptm-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.12);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.iptm-bg {
    background: linear-gradient(135deg, #ff9a9e 0%, #fecfef 100%);
    color: #333;
}

.metrics-tile-back.iptm-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.15);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.iptm-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.12);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.sa-bg {
    background: linear-gradient(135deg, #a18cd1 0%, #fbc2eb 100%);
    color: #fff;
}

.metrics-tile-back.sa-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.sa-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.admet-bg {
    background: linear-gradient(135deg, #fad0c4 0%, #ffd1ff 100%);
    color: #333;
}

.metrics-tile-back.admet-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.15);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.admet-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.12);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.logp-bg {
    background: linear-gradient(135deg, #ffecd2 0%, #fcb69f 100%);
    color: #333;
}

.metrics-tile-back.logp-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.15);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.logp-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.12);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.solubility-bg {
    background: linear-gradient(135deg, #a8edea 0%, #fed6e3 100%);
    color: #333;
}

.metrics-tile-back.solubility-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.15);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.solubility-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.12);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.toxicity-bg {
    background: linear-gradient(135deg, #ff9a9e 0%, #fecfef 100%);
    color: #333;
}

.metrics-tile-back.toxicity-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.15);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.toxicity-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.12);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.targets-bg {
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
    color: #fff;
}

.metrics-tile-back.targets-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.targets-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.weight-bg {
    background: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
    color: #fff;
}

.metrics-tile-back.weight-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.weight-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.druglikeness-bg {
    background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
    color: #fff;
}

.metrics-tile-back.druglikeness-bg::before {
    content: "";
    position: absolute;
    top: -40px;
    right: -40px;
    width: 80px;
    height: 80px;
    background: rgba(255, 255, 255, 0.1);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-back.druglikeness-bg::after {
    content: "";
    position: absolute;
    bottom: -20px;
    left: -20px;
    width: 40px;
    height: 40px;
    background: rgba(255, 255, 255, 0.08);
    border-radius: 50%;
    z-index: -1;
}

.metrics-tile-container * {
    transition: all 0.3s ease;
}

@keyframes pulse-glow {
    0% { box-shadow: 0 6px 20px rgba(0,0,0,0.1); }
    50% { box-shadow: 0 6px 25px rgba(0,0,0,0.2); }
    100% { box-shadow: 0 6px 20px rgba(0,0,0,0.1); }
}

@keyframes tile-entrance {
    from {
        opacity: 0;
        transform: translateY(20px) scale(0.95);
    }
    to {
        opacity: 1;
        transform: translateY(0) scale(1);
    }
}

@media (max-width: 900px) {
    .metrics-tile-container {
        height: 140px;
    }
    
    .metrics-tile-back {
        font-size: 13px;
        padding: 12px;
        line-height: 0.5;
    }
    
    .metrics-tile-back h4 {
        font-size: 16px;
        margin: 0 0 6px 0;
        line-height: 0.5;
    }
    
    .metrics-tile-back div {
        font-size: 13px;
        line-height: 0.5;
    }
}
</style>
""" 