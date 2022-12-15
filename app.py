import vwtds
from flask import Flask


app = Flask(__name__)

@app.route('/', methods=['GET'])
def get_vwtds():
    results = vwtds.run_vwtds()
    return results

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=2000)
