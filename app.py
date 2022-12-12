import vwtds
from flask import Flask


app = Flask(__name__)

@app.route('/get_vwtds_data', methods=['GET'])
def get_vwtds():
    results = vwtds.run_vwtds()
    return results

if __name__ == '__main__':
    app.run()