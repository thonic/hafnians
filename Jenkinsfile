pipeline {
    agent any
    stages {
        stage('Run hafnian.py') {
            when { allOf {
                branch 'master';
                changelog '^RUN .*'
                }
            }
            steps {
                sh '''
                    pwd;
                    m=$(git log -1 --pretty=%B | sed 's/RUN \\(.*\\)/\\1/')
                    echo $m
                    # /home/thonic/anaconda3/envs/partial-derivatives/bin/python /home/thonic/git/partial-derivatives/hafnian.py
                '''
            }
        }
    }
}