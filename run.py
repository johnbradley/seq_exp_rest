import os
import seq_exp.seq_exp as seq_exp

if __name__ == "__main__":
    url = os.environ['DBURL']
    if not url:
        url='sqlite:///:memory:'
    app, db = seq_exp.setup_api_and_db(url)
    app.run(debug=True)
