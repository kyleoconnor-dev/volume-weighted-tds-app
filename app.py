import vwtds
from flask import Flask
from flask_cors import CORS


app = Flask(__name__)
CORS(app)

@app.route('/', methods=['GET'])
def get_vwtds():
    results = vwtds.run_vwtds()
    return results

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=2000)
